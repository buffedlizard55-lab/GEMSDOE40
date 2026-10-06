#!/usr/bin/env python3
"""Generate the preregistered H8 ASA/SPI depth-cluster candidate (frozen settings).

Reads only: the pinned feature raster (bands by name), the sample template
(shape/CRS/transform/footprint), the publicly supplied known-fault raster
(exact-pixel mask + flank ramp), and this repo's committed H4 Euler solution
CSV as corroborating evidence. Never reads proxy truth, prior predictions,
holdout labels, or scores.
"""
from __future__ import annotations

import argparse
import csv
import datetime as _dt
import gzip
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe40 import asa_spi  # noqa: E402
from gemsdoe40.raster import (  # noqa: E402
    band_index_by_name,
    canonical_pixel_sha256,
    write_candidate,
    write_json,
)

# --------------------------- frozen configuration ---------------------------
CONFIG = {
    "id": "H8",
    "pad_cells": 128,
    "peak_quantile": 0.95,
    "magnetic": {
        "band_name": "tmi",
        "continuation_m": 300.0,
        "margin_m": 1500.0,
        "depth_bounds_m": [100.0, 5000.0],
    },
    "gravity": {
        "band_name": "iso_grav_anom",
        "continuation_m": 500.0,
        "margin_m": 3000.0,
        "depth_bounds_m": [100.0, 8000.0],
    },
    "corroboration_xy_m": 300.0,
    "corroboration_depth_floor_m": 300.0,
    "corroboration_depth_rel": 0.4,
    "weights": {
        "depth_scale_m": 1500.0,
        "std_floor_m": 250.0,
        "std_rel": 0.35,
        "corroborated_boost": 2.5,
        "non_persistent_factor": 0.35,
    },
    "kde": {"sigma_px": 2.0, "truncate": 3.0},
    "support_floor": 0.02,
    "catalogue_ramp_m": [150.0, 450.0],
    "corroboration_csv": "docs/downloads/h4-euler-solutions.csv.gz",
    "preregistration": "docs/research/h8-preregistration-20261006.md",
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_family(name: str, field: np.ndarray, valid: np.ndarray, h4_solutions: dict) -> tuple[np.ndarray, dict]:
    """Registered dense construction (preregistration item-7 amendment).

    Every in-footprint cell with a valid contact depth and analytic-signal
    amplitude above the registered floor is a weighted solution; a Gaussian
    KDE smooths the weighted solution density into the continuous field.
    """
    fcfg = CONFIG[name]
    t0 = time.time()
    stats: dict = {"family": name}
    base = asa_spi.analytic_signal_depth(field, pad_cells=CONFIG["pad_cells"])
    cont_field = asa_spi.upward_continue(field, fcfg["continuation_m"], pad_cells=CONFIG["pad_cells"])
    cont = asa_spi.analytic_signal_depth(cont_field, pad_cells=CONFIG["pad_cells"])
    del cont_field
    a, z = base["a_smooth"], base["z"]
    zc = cont["z"] - fcfg["continuation_m"]
    far = valid & asa_spi.far_from_invalid(valid, fcfg["margin_m"])
    zmin, zmax = fcfg["depth_bounds_m"]
    depth_ok = np.isfinite(z) & (z >= zmin) & (z <= zmax)
    base_set = far & depth_ok
    thresh = float(np.quantile(a[base_set], CONFIG["peak_quantile"])) if base_set.any() else np.inf
    sol = base_set & (a >= thresh)
    stats["solution_cells"] = int(sol.sum())
    stats["amplitude_threshold"] = thresh
    stats["depth_in_bounds_cells"] = int(depth_ok.sum())

    # Euler corroboration on solution cells (same family)
    sol_key = "tmi" if name == "magnetic" else "gravity"
    hh = h4_solutions[sol_key]
    corroborated = np.zeros(field.shape, dtype=bool)
    if sol.any() and hh["row"].size:
        srow, scol = np.nonzero(sol)
        from scipy.spatial import cKDTree
        tree = cKDTree(np.column_stack([hh["row"], hh["col"]]) * asa_spi.RESOLUTION_M)
        dist, idx = tree.query(np.column_stack([srow, scol]) * asa_spi.RESOLUTION_M, k=1,
                               distance_upper_bound=CONFIG["corroboration_xy_m"])
        found = np.isfinite(dist)
        if found.any():
            tol = np.maximum(CONFIG["corroboration_depth_floor_m"],
                             CONFIG["corroboration_depth_rel"] * np.abs(z[srow, scol]))
            cand = np.where(found, idx, 0)
            ok = found & (np.abs(hh["depth_m"][cand] - z[srow, scol]) <= tol)
            corroborated[srow[ok], scol[ok]] = True
    stats["cells_corroborated_by_h4"] = int(corroborated.sum())

    weight, persist = asa_spi.dense_family_weights(
        z, sol, cont["a_smooth"], zc, corroborated,
        depth_scale_m=CONFIG["weights"]["depth_scale_m"],
        std_floor_m=CONFIG["weights"]["std_floor_m"],
        std_rel=CONFIG["weights"]["std_rel"],
        non_persistent_factor=CONFIG["weights"]["non_persistent_factor"],
        corroborated_boost=CONFIG["weights"]["corroborated_boost"],
    )
    stats["cells_persistent"] = int(persist.sum())
    stats["cells_weighted_nonzero"] = int((weight > 0).sum())
    from scipy import ndimage as _ndi
    kde = _ndi.gaussian_filter(weight, sigma=CONFIG["kde"]["sigma_px"],
                               truncate=CONFIG["kde"]["truncate"], mode="constant", cval=0.0)
    stats["seconds"] = round(time.time() - t0, 1)

    row, col = np.nonzero(weight > 0)
    cloud = {
        "family": name, "row": row.astype(np.float64), "col": col.astype(np.float64),
        "depth_m": z[row, col], "amplitude": a[row, col],
        "persistent": persist[row, col],
        "corroborated": corroborated[row, col], "weight": weight[row, col],
    }
    return kde, {"stats": stats, "cloud": cloud}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, default=ROOT / "work/h8")
    parser.add_argument("--out-dir", type=Path, default=ROOT / "docs/downloads",
                        help="artifact directory (reproduction runs use work/h8/repro)")
    parser.add_argument("--stamp", type=str, default=None,
                        help="UTC timestamp token for the artifact name; default now()")
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    data = ROOT / "data"
    features, sample, labels = (
        data / "training_features.tif", data / "sample_submission.tif", data / "labels.tif")
    pins = {
        "training_features.tif": "4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5",
        "sample_submission.tif": "2176d08e485aa2cd2860ce8df539db4faf4d76163b38a4dd8c30a40454d35cbc",
        "labels.tif": "7ba308ccdc4418b31a178f4f1ef21aaa6e152e4028f2f6f64b01f7eb25ae4093",
    }
    input_hashes = {}
    for rel, expected in pins.items():
        got = sha256_file(ROOT / "data" / rel)
        if got != expected:
            raise RuntimeError(f"pinned input changed: {rel} {got}")
        input_hashes[rel] = got

    stamp = args.stamp or _dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    with rasterio.open(features) as ds:
        tmi = ds.read(band_index_by_name(ds, CONFIG["magnetic"]["band_name"])).astype(np.float32)
        grav = ds.read(band_index_by_name(ds, CONFIG["gravity"]["band_name"])).astype(np.float32)
    with rasterio.open(sample) as ds:
        footprint = np.isfinite(ds.read(1))
        transform = ds.transform
        crs = ds.crs
    with rasterio.open(labels) as ds:
        labels_arr = ds.read(1)
    # NB: the feature raster uses a finite -3.4e38 sentinel; plain isfinite() is wrong here.
    valid = footprint & asa_spi._finite_mask(tmi) & asa_spi._finite_mask(grav)
    if not valid.any():
        raise RuntimeError("no valid cells in the footprint")

    h4_solutions = asa_spi.load_h4_solutions(ROOT / CONFIG["corroboration_csv"],
                                             transform_origin=(transform.c, transform.f))

    kde_mag, mag = build_family("magnetic", tmi, valid, h4_solutions)
    del tmi
    # Registered gravity family: the pipeline runs on Gz (first vertical derivative
    # of the isostatic anomaly) as a local top-edge approximation, and the
    # continuation height applies to that Gz field.
    gz = asa_spi.vertical_derivative(grav, pad_cells=CONFIG["pad_cells"])
    del grav
    kde_grav, gravf = build_family("gravity", gz, valid, h4_solutions)
    del gz

    fused = asa_spi.fuse_families(kde_mag, kde_grav)
    pos = fused[fused > 0]
    if pos.size == 0:
        raise RuntimeError("H8 field is empty; refusing to publish")
    field = np.clip(fused / np.quantile(pos, 0.995), 0.0, 1.0)
    ramp = asa_spi.catalogue_ramp(labels_arr, zero_m=CONFIG["catalogue_ramp_m"][0],
                                  full_m=CONFIG["catalogue_ramp_m"][1])
    field = field * ramp
    field[(labels_arr == 1) & footprint] = 0.0
    field = np.where(field < CONFIG["support_floor"], 0.0, field)
    if (field > 0).sum() == 0:
        raise RuntimeError("H8 support vanished after masking; refusing to publish")
    field = field / field.max()
    field = np.where(footprint, field, 0.0)

    cloud_hash_source = np.concatenate([
        mag["cloud"]["row"], mag["cloud"]["col"], mag["cloud"]["depth_m"],
        gravf["cloud"]["row"], gravf["cloud"]["col"], gravf["cloud"]["depth_m"]])
    cloud_digest = hashlib.sha256(np.ascontiguousarray(cloud_hash_source.astype("<f8")).tobytes()).hexdigest()[:12]
    basename = f"gemsdoe40-h8-asa-spi-depthkde-{stamp}-{cloud_digest}"

    out_nan = out_dir / f"{basename}-nan.tif"
    write_candidate(out_nan, field.astype(np.float32), sample,
                    description="GEMSDOE40 preregistered H8 analytic-signal local-wavenumber depth-cluster KDE")
    zeros = field.astype(np.float32).copy()
    zeros[~footprint] = 0.0
    out_zeros = out_dir / f"{basename}-zeros.tif"
    with rasterio.open(sample) as ds:
        profile = ds.profile.copy()
    profile.update(driver="GTiff", count=1, dtype="float32", compress="DEFLATE", predictor=3,
                   BIGTIFF="IF_SAFER")
    profile.pop("nodata", None)
    with rasterio.open(out_zeros, "w", **profile) as dst:
        dst.write(zeros, 1)
        dst.update_tags(1, description="GEMSDOE40 H8 portal-safe twin: all finite, zeros outside footprint")
    import zipfile
    out_zip = out_dir / f"{basename}-zeros.zip"
    zip_time = _dt.datetime.strptime(stamp, "%Y%m%dT%H%M%SZ").timetuple()[:6]
    with zipfile.ZipFile(out_zip, "w", zipfile.ZIP_DEFLATED) as archive:
        info = zipfile.ZipInfo(out_zeros.name, date_time=zip_time)
        info.compress_type = zipfile.ZIP_DEFLATED
        with open(out_zeros, "rb") as fh:
            archive.writestr(info, fh.read())

    # compressed point cloud
    cloud_csv = out_dir / f"{basename}-cloud.csv.gz"
    import io
    with open(cloud_csv, "wb") as raw, io.TextIOWrapper(gzip.GzipFile(
            filename="", mode="wb", fileobj=raw, mtime=0), newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["family", "easting_m", "northing_m", "effective_depth_m",
                         "amplitude", "persistent", "corroborated", "weight"])
        for fam, info in (("magnetic", mag), ("gravity", gravf)):
            c = info["cloud"]
            for i in range(c["row"].size):
                easting = transform.c + c["col"][i] * transform.a
                northing = transform.f + c["row"][i] * transform.e
                writer.writerow([fam, f"{easting:.2f}", f"{northing:.2f}", f"{c['depth_m'][i]:.3f}",
                                 f"{c['amplitude'][i]:.6g}", int(c["persistent"][i]),
                                 int(c["corroborated"][i]), f"{c['weight'][i]:.6g}"])

    code_hashes = {rel: sha256_file(ROOT / rel) for rel in (
        "src/gemsdoe40/asa_spi.py", "scripts/run_h8_asa_spi.py")}
    receipt = {
        "basename": basename,
        "stamp": stamp,
        "cloud_digest": cloud_digest,
        "config": CONFIG,
        "input_sha256": input_hashes,
        "code_sha256": code_hashes,
        "preregistration_sha256": sha256_file(ROOT / CONFIG["preregistration"]),
        "magnetic": mag["stats"],
        "gravity": gravf["stats"],
        "output": {
            "nan_tif": {"path": str(out_nan.relative_to(ROOT)) if out_nan.is_relative_to(ROOT) else str(out_nan), "sha256": sha256_file(out_nan),
                        "bytes": out_nan.stat().st_size},
            "zeros_tif": {"path": str(out_zeros.relative_to(ROOT)), "sha256": sha256_file(out_zeros),
                          "bytes": out_zeros.stat().st_size},
            "zip": {"path": str(out_zip.relative_to(ROOT)), "sha256": sha256_file(out_zip),
                    "bytes": out_zip.stat().st_size},
            "cloud_csv": {"path": str(cloud_csv.relative_to(ROOT)), "sha256": sha256_file(cloud_csv),
                          "bytes": cloud_csv.stat().st_size},
        },
        "field_stats": {
            "positive_cells": int((field > 0).sum()),
            "distinct_positive_values": int(np.unique(field[field > 0]).size),
            "confidence_mass_sum": float(field.sum()),
            "max": float(field.max()),
        },
    }
    write_json(work / "generation.json", receipt)
    np.savez_compressed(work / "field.npz", field=field.astype(np.float32), footprint=footprint)
    print(f"PASS H8 generation: {basename}")
    print(json.dumps(receipt["field_stats"], indent=2))
    print(json.dumps({"magnetic": mag["stats"], "gravity": gravf["stats"]}, indent=2))


if __name__ == "__main__":
    main()
