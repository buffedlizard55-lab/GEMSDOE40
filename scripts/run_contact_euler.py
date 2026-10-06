#!/usr/bin/env python3
"""Generate the locked H4 cloud/KDE, without reading holdout or prior predictions.

Writes to ignored work/ first. Publication occurs only after the independent
format + novelty audit; no competition upload code exists in this command.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import gc
import gzip
from importlib.metadata import version
import json
from pathlib import Path
import sys
import time

import numpy as np
import rasterio
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gemsdoe40.contact_euler import (  # noqa: E402
    FIELDS, SETTINGS, config_dict, spectral_gradients, solution_cloud,
    cluster_weights, family_kde, combine_families,
)
from gemsdoe40.raster import band_index_by_name, validate_candidate, write_candidate, write_json  # noqa: E402
from acquire_data import PINS, digest  # noqa: E402


def save_cloud(path: Path, clouds: dict, transform) -> None:
    """Keep all QC-passing solutions, including cluster weight=0, for audit."""
    with path.open("wb") as stream:
        with gzip.GzipFile(filename="", fileobj=stream, mode="wb", mtime=0) as zipped:
            import io
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                out = csv.writer(text)
                fields = ["family", "window_cells", "easting_m", "northing_m", "effective_depth_m",
                          "conditional_depth_se_m", "relative_residual", "design_condition", "contact_offset",
                          "identifiable_rank", "cluster_weight", "matching_neighbours", "depth_mad_m"]
                out.writerow(fields)
                for name, cloud in clouds.items():
                    for c in cloud:
                        east, north = transform * (float(c["col"]) + .5, float(c["row"]) + .5)
                        out.writerow([name, int(c["window"]), f"{east:.4f}", f"{north:.4f}",
                                      *[f"{float(c[f]):.9g}" for f in ("depth_m", "depth_se_m", "residual", "condition", "offset")],
                                      int(c["rank"]), f"{c['weight']:.9g}", int(c["neighbours"]), f"{c['depth_mad_m']:.9g}"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data")
    parser.add_argument("--output", type=Path, default=ROOT / "work/h4")
    args = parser.parse_args()
    args.data = args.data.resolve()
    args.output = args.output.resolve()
    if not args.output.is_relative_to(ROOT / "work"):
        parser.error("generation output must stay under the repository's ignored work/ directory")
    started = time.monotonic()
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = {"experiment": "H4 offset-aware rank-adaptive contact Euler depth-clustering",
                "started_utc": datetime.now(timezone.utc).isoformat(), "configuration": config_dict(),
                "software": {"python": sys.version.split()[0], **{p: version(p) for p in ("numpy", "scipy", "rasterio", "affine")}},
                "compute": "CPU; OPENBLAS_NUM_THREADS=1; FFT workers=2; no learned weights or GPU",
                "construction_uses_proxy_or_prior_predictions": False, "input_files": {},
                "preregistration_sha256": digest(ROOT / "docs/research/h4-preregistration-20261006.md"),
                "code_sha256": {str(p.relative_to(ROOT)): digest(p) for p in (Path(__file__), ROOT / "src/gemsdoe40/contact_euler.py")}}
    for name in ("training_features.tif", "sample_submission.tif", "labels.tif"):
        p = args.data / name
        sha = digest(p)
        if sha != PINS[name]:
            raise RuntimeError(f"input hash mismatch: {name}")
        manifest["input_files"][name] = {"sha256": sha, "bytes": p.stat().st_size}
    sample = args.data / "sample_submission.tif"
    with rasterio.open(sample) as ds:
        footprint = np.isfinite(ds.read(1))
        shape, crs, transform = ds.shape, ds.crs, ds.transform
    with rasterio.open(args.data / "labels.tif") as ds:
        if ds.shape != shape or ds.crs != crs or ds.transform != transform:
            raise RuntimeError("labels do not match sample grid")
        known = (ds.read(1) == 1) & footprint
    manifest["shape"] = list(shape)
    manifest["footprint_pixels"] = int(footprint.sum())
    manifest["exact_known_mask_pixels"] = int(known.sum())
    manifests, clouds, kdes = {}, {}, {}
    for config in FIELDS:
        t0 = time.monotonic()
        with rasterio.open(args.data / "training_features.tif") as ds:
            if ds.shape != shape or ds.crs != crs or ds.transform != transform:
                raise RuntimeError("features do not match sample grid")
            band = band_index_by_name(ds, config.name)
            field = ds.read(band)
            valid = (ds.read_masks(band) > 0) & footprint & np.isfinite(field) & (np.abs(field) < 1e30)
        print(f"{config.name}: band {band}, {int(valid.sum()):,} real cells; spectral derivatives", flush=True)
        tx, ty, tz = spectral_gradients(field, valid, cell_m=100, height_m=config.continuation_m,
                                        derivative_order=config.derivative_order, pad=SETTINGS["fft_pad"])
        # Padding with invalid cells makes the rectangular raster border explicit.
        distance = ndimage.distance_transform_edt(np.pad(valid, 1, constant_values=False))[1:-1, 1:-1] * 100.0
        safe = valid & (distance >= config.edge_guard_m)
        del distance, field
        field_clouds, receipts = [], []
        for win in config.windows:
            print(f"{config.name}: solving {win}x{win} windows", flush=True)
            cloud, receipt = solution_cloud(tx, ty, tz, safe, config, win, stride=SETTINGS["stride"])
            field_clouds.append(cloud); receipts.append(receipt)
            print(json.dumps(receipt), flush=True)
        del tx, ty, tz
        gc.collect()
        cloud = cluster_weights(np.concatenate(field_clouds))
        kde, summary = family_kde(cloud, footprint)
        clouds[config.name] = cloud
        kdes[config.name] = kde
        np.save(args.output / f"{config.name}-kde.npy", kde)
        np.savez_compressed(args.output / f"{config.name}-solutions.npz", cloud=cloud)
        manifests[config.name] = {"band_index": band, "real_input_cells": int(valid.sum()),
                                  "guarded_cells": int(safe.sum()), "windows": receipts, "kde": summary,
                                  "elapsed_s": round(time.monotonic() - t0, 3)}
        print(f"{config.name} clustered: {json.dumps(summary)}", flush=True)
        write_json(args.output / "partial-generation.json", {**manifest, "families": manifests})
    prediction = combine_families(kdes["tmi"], kdes["iso_grav_anom"], footprint, known)
    temp = args.output / "h4.tif"
    write_candidate(temp, prediction, sample, description="GEMSDOE40 H4 SI0+A rank-aware TMI + dG/dz depth-consistent continuous Euler KDE; research only")
    receipt = validate_candidate(temp, sample)
    if not receipt["valid"]:
        raise RuntimeError(f"format audit failed: {receipt['issues']}")
    filename = f"gemsdoe40-h4-contact-offset-depthkde-20261006-{receipt['canonical_pixels_sha256'][:12]}.tif"
    destination = args.output / filename
    if destination.exists() and digest(destination) != receipt["sha256"]:
        raise RuntimeError("refuse overwrite of different bytes with same candidate name")
    temp.replace(destination)
    receipt["path"] = str(destination.relative_to(ROOT))
    cloud_path = args.output / "h4-euler-solutions.csv.gz"
    save_cloud(cloud_path, clouds, transform)
    manifest.update(families=manifests, candidate_format=receipt,
                    cloud={"path": str(cloud_path.relative_to(ROOT)), "sha256": digest(cloud_path),
                           "bytes": cloud_path.stat().st_size, "records": sum(len(c) for c in clouds.values())},
                    positive_pixels=int(np.count_nonzero(prediction[footprint])),
                    unique_finite_prediction_values=int(np.unique(prediction[footprint]).size),
                    in_footprint_mass=float(prediction[footprint].sum(dtype=np.float64)),
                    completed_utc=datetime.now(timezone.utc).isoformat(),
                    elapsed_s=round(time.monotonic() - started, 3),
                    state="GENERATED; NOT YET NOVELTY/PROXY AUDITED; NO UPLOAD")
    write_json(args.output / "generation.json", manifest)
    print(json.dumps({"candidate": receipt["path"], "format_pass": receipt["valid"],
                      "positive_cells": manifest["positive_pixels"], "elapsed_s": manifest["elapsed_s"]}), flush=True)


if __name__ == "__main__":
    main()
