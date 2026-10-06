#!/usr/bin/env python3
"""Generate the H8 candidate into ignored ``work/`` only.

Reads nothing but the hash-pinned competition rasters.  No prior prediction and
no holdout/label raster enters generation; the labels are read by the separate
audit/validation runners, never here.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from importlib.metadata import version
import json
from pathlib import Path
import sys
import time

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gemsdoe40.contact_euler import spectral_gradients  # noqa: E402
from gemsdoe40.emission import value_ranked_thinning  # noqa: E402
from gemsdoe40.h8_euler import (FAMILIES, SETTINGS, combine_weights, config_dict,  # noqa: E402
                                cross_family_corroboration, depth_consensus_weights,
                                kde_field, lineament_coherence)
from gemsdoe40.h8_solver import solve_window  # noqa: E402
from gemsdoe40.raster import validate_candidate, write_candidate, write_json  # noqa: E402
from acquire_data import PINS, digest  # noqa: E402

SENTINEL = -3.4028234663852886e38


def read_band(path: Path, name: str):
    with rasterio.open(path) as src:
        index = None
        for i in range(1, src.count + 1):
            tag = (src.descriptions[i - 1] or "").split(" - ")[0].strip()
            if tag == name:
                index = i
                break
        if index is None:
            raise RuntimeError(f"band {name!r} not found in {path}")
        band = src.read(index).astype(np.float32)
    valid = np.isfinite(band) & (band > SENTINEL / 2.0)
    band = np.where(valid, band, np.nan)
    return band, valid


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data")
    parser.add_argument("--output", type=Path, default=ROOT / "work" / "h8")
    args = parser.parse_args()
    args.output = args.output.resolve()
    if not args.output.is_relative_to(ROOT / "work"):
        parser.error("generation output must stay under ignored work/")
    args.output.mkdir(parents=True, exist_ok=True)

    started = time.monotonic()
    manifest = {
        "experiment": "H8 lineament-weighted cross-family-corroborated SI=0 Euler depth-clustering",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "configuration": config_dict(),
        "software": {"python": sys.version.split()[0],
                     **{p: version(p) for p in ("numpy", "scipy", "rasterio")}},
        "compute": "CPU; no GPU, no learned weights",
        "construction_reads_labels_or_priors": False,
        "input_files": {},
    }
    for name in ("training_features.tif", "sample_submission.tif", "labels.tif"):
        path = args.data / name
        sha = digest(path)
        if path.name in PINS and sha != PINS[path.name]:
            raise RuntimeError(f"input hash mismatch: {name}")
        manifest["input_files"][name] = {"sha256": sha, "bytes": path.stat().st_size}

    with rasterio.open(args.data / "sample_submission.tif") as ds:
        footprint = np.isfinite(ds.read(1))
        shape, crs, transform = ds.shape, ds.crs, ds.transform
    manifest["grid"] = {"shape": list(shape), "crs": crs.to_string(),
                        "transform": list(transform)[:6], "footprint_cells": int(footprint.sum())}

    clouds = {}
    totals = []
    for family_id, family in enumerate(FAMILIES):
        band, valid = read_band(args.data / "training_features.tif", family.name)
        manifest.setdefault("bands", {})[family.name] = {
            "family_id": family_id, "valid_cells": int(valid.sum()),
            "continuation_m": family.continuation_m, "derivative_order": family.derivative_order,
            "windows": list(family.windows), "max_depth_m": family.max_depth_m}
        pieces = []
        for window in family.windows:
            tx, ty, tz = spectral_gradients(band, valid, height_m=family.continuation_m,
                                            derivative_order=family.derivative_order,
                                            pad=SETTINGS["fft_pad"])
            cloud, stats = solve_window(tx, ty, tz, valid, family, window,
                                        family_id=family_id, stride=SETTINGS["stride"])
            totals.append(stats)
            pieces.append(cloud)
            print(f"  {family.name} w={window}: accepted {stats['accepted']} "
                  f"of {stats['windows_tested']} windows", flush=True)
            del tx, ty, tz
        clouds[family.name] = np.concatenate(pieces) if pieces else np.empty(0, dtype=clouds_dtype_placeholder())
        del band, valid

    mag = clouds["rtp"]
    grav = clouds["iso_grav_anom"]
    for cloud in (mag, grav):
        if not len(cloud):
            raise RuntimeError("empty solution cloud; refusing to continue")
    for cloud in (mag, grav):
        cloud[:] = depth_consensus_weights(cloud, xy_radius_m=SETTINGS["cluster_xy_m"],
                                           depth_floor_m=SETTINGS["cluster_depth_floor_m"],
                                           depth_fraction=SETTINGS["cluster_depth_fraction"],
                                           min_neighbours=SETTINGS["min_other_neighbours"],
                                           depth_decay_m=SETTINGS["depth_decay_m"])
        cloud["coherence"] = lineament_coherence(cloud, xy_radius_m=SETTINGS["coherence_xy_m"],
                                                min_neighbours=SETTINGS["min_coherence_neighbours"])
    mag_flag, grav_flag = cross_family_corroboration(mag, grav,
                                                     xy_radius_m=SETTINGS["corroboration_xy_m"],
                                                     depth_m=SETTINGS["corroboration_depth_m"])
    mag["corroborated"], grav["corroborated"] = mag_flag, grav_flag
    weights = {"rtp": combine_weights(mag), "iso_grav_anom": combine_weights(grav)}

    # one joint cloud: every retained solution from both families, in one KDE
    joint = np.concatenate([mag, grav])
    joint_weight = np.concatenate([weights["rtp"], weights["iso_grav_anom"]])
    field, kde_stats = kde_field(joint, joint_weight, footprint)

    emitted = value_ranked_thinning(field, footprint, spacing_px=SETTINGS["emission_spacing_px"],
                                    max_mass=SETTINGS["mass_budget"], min_score=0.0)
    support_mask = emitted > 0
    prediction = np.where(support_mask, field, 0.0).astype(np.float32)
    prediction[~footprint] = 0.0
    if not np.isfinite(prediction[footprint]).all():
        raise RuntimeError("non-finite prediction inside footprint")
    mass = int(support_mask.sum())
    if mass == 0:
        raise RuntimeError("empty emission")

    candidate = args.output / "h8-euler-depthcluster-candidate.tif"
    write_candidate(candidate, prediction, args.data / "sample_submission.tif",
                    description="GEMSDOE40 H8 lineament-weighted cross-family SI=0 Euler depth-cluster KDE")
    receipt = validate_candidate(candidate, args.data / "sample_submission.tif")

    per_family = {}
    for name in clouds:
        cloud, weight = clouds[name], weights[name]
        retained = weight > 0
        per_family[name] = {
            "accepted_solutions": int(len(cloud)),
            "clustered_solutions": int(retained.sum()),
            "weight_sum": float(weight.sum()),
            "depth_m_p05_p50_p95": (np.percentile(cloud["depth_m"][retained], [5, 50, 95]).tolist()
                                    if retained.any() else None),
            "corroborated_fraction": float((cloud["corroborated"][retained] > 0).mean()) if retained.any() else None,
            "coherence_p50": float(np.median(cloud["coherence"][retained])) if retained.any() else None,
            "median_neighbours": float(np.median(cloud["neighbours"][retained])) if retained.any() else None,
        }
    values = prediction[support_mask].astype(np.float64)
    document = {
        "experiment": manifest["experiment"],
        "started_utc": manifest["started_utc"],
        "elapsed_s": round(time.monotonic() - started, 1),
        "configuration": manifest["configuration"],
        "software": manifest["software"],
        "input_files": manifest["input_files"],
        "bands": manifest.get("bands", {}),
        "grid": manifest["grid"],
        "window_totals": totals,
        "families": per_family,
        "kde": kde_stats,
        "emission": {
            "spacing_px": SETTINGS["emission_spacing_px"],
            "mass_budget": SETTINGS["mass_budget"],
            "emitted_pixels": mass,
            "distinct_values": int(np.unique(values).size),
            "value_min": float(values.min()), "value_max": float(values.max()),
            "value_p50": float(np.median(values)),
            "value_fraction_at_one": float((values >= 0.999999).mean()),
        },
        "candidate": receipt,
        "construction_reads_labels_or_priors": False,
        "not_a_score": "Every field here is a measurement or a model; no organizer score exists for this file.",
    }
    write_json(args.output / "h8-generation.json", document)
    np.save(args.output / "h8-field.npy", field)
    print(f"PASS: {mass} emitted pixels; {receipt['candidate']['in_footprint_min']}.."
          f"{receipt['candidate']['in_footprint_max']}; {candidate}")
    print(f"  sha256 {receipt['sha256']}")
    print(f"  pixel sha256 {receipt['canonical_pixels_sha256']}")
    return 0


def clouds_dtype_placeholder():
    from gemsdoe40.h8_euler import CLOUD_DTYPE
    return CLOUD_DTYPE


if __name__ == "__main__":
    sys.exit(main())
