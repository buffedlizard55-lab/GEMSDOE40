#!/usr/bin/env python3
"""Build the locked H4 potential-field Euler/3-D-cluster KDE candidate.

This construction stage reads only the hash-pinned feature stack, sample grid, and exact known-label
mask (for exact-pixel suppression). It deliberately does not open the proxy, any holdout truth,
prior prediction, or leaderboard data. Holdout evaluation is a separate command.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import rasterio

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gems40.pins import PINNED_FILES
from gemsdoe40.euler_h4 import (
    build_h4_prediction,
    potential_field_derivatives,
    solve_euler_cloud,
)
from gemsdoe40.raster import band_index_by_name, canonical_pixel_sha256, validate_candidate, write_candidate, write_json


INPUT_SHA256 = {
    "features": PINNED_FILES["training_features.tif"]["sha256"],
    "sample": PINNED_FILES["example_submission.tif"]["sha256"],
    "labels": PINNED_FILES["existing_faults.tif"]["sha256"],
}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_inputs(paths: dict[str, Path]) -> dict[str, str]:
    verified: dict[str, str] = {}
    for key, path in paths.items():
        actual = file_sha256(path)
        expected = INPUT_SHA256[key]
        if actual != expected:
            raise ValueError(f"{key} SHA-256 mismatch: {actual} != pinned {expected}")
        verified[key] = actual
    return verified


def _assert_aligned(feature_ds, sample_ds, labels_ds) -> None:
    reference = (sample_ds.shape, sample_ds.crs, sample_ds.transform)
    for ds, name in ((feature_ds, "features"), (labels_ds, "labels")):
        if (ds.shape, ds.crs, ds.transform) != reference:
            raise ValueError(f"{name} raster is not exactly aligned to sample_submission.tif")
    if feature_ds.count != 19 or labels_ds.count != 1 or sample_ds.count != 1:
        raise ValueError("expected 19 feature bands and single-band sample/labels rasters")
    if not sample_ds.crs or sample_ds.crs.to_epsg() != 32611:
        raise ValueError(f"sample CRS must be EPSG:32611, got {sample_ds.crs}")
    if abs(sample_ds.transform.a - 100.0) > 1e-9 or abs(sample_ds.transform.e + 100.0) > 1e-9:
        raise ValueError("sample grid must be 100 m north-up")


def _read_band(feature_ds, name: str) -> tuple[np.ndarray, np.ndarray, int]:
    band = band_index_by_name(feature_ds, name)
    data = feature_ds.read(band).astype(np.float32, copy=False)
    valid = np.isfinite(data) & (np.abs(data) < 1e30)
    return data, valid, band


def build(args: argparse.Namespace) -> tuple[Path, dict]:
    paths = {"features": args.features, "sample": args.sample, "labels": args.labels}
    verified = verify_inputs(paths)
    with rasterio.open(args.features) as features, rasterio.open(args.sample) as sample, rasterio.open(args.labels) as labels_ds:
        _assert_aligned(features, sample, labels_ds)
        template = sample.read(1)
        footprint = np.isfinite(template)
        labels = labels_ds.read(1)
        if not np.isin(np.unique(labels), [-1, 0, 1]).all():
            raise ValueError("known-label raster has unexpected values")
        known = (labels == 1) & footprint
        transform = features.transform
        bands = {
            name: band_index_by_name(features, name)
            for name in ("rtp", "iso_grav_anom")
        }
        clouds = {}
        solve_stats = {}
        field_stats = {}
        for name, si in (("rtp", 0.0), ("iso_grav_anom", -1.0)):
            raw, band_valid, band = _read_band(features, name)
            valid = footprint & band_valid
            if not valid.any():
                raise RuntimeError(f"no valid in-footprint pixels for {name}")
            filled, derivative_valid, tx, ty, tz = potential_field_derivatives(raw)
            derivative_valid &= valid
            field_stats[name] = {
                "band": int(band),
                "valid_cells_in_sample_footprint": int(valid.sum()),
                "field_min_max": [float(raw[valid].min()), float(raw[valid].max())],
                "field_std": float(raw[valid].std(dtype=np.float64)),
                "derivative_convention": "100 m central differences east/north; downward-positive +|k| vertical derivative, reflect-padded by 128 cells",
            }
            cloud, counts = solve_euler_cloud(
                filled,
                tx,
                ty,
                tz,
                derivative_valid,
                transform=transform,
                si=si,
                field=name,
            )
            if len(cloud) == 0:
                raise RuntimeError(
                    f"H4 stop: the locked {name} Euler solve produced no accepted solutions; "
                    f"window-gate counts={counts}; no threshold relaxation or gradient fallback is permitted"
                )
            clouds[name] = cloud
            solve_stats[name] = counts
            # Release full-grid arrays before processing the second field.
            del raw, valid, filled, derivative_valid, tx, ty, tz

        prediction, h4_summary = build_h4_prediction(
            clouds["rtp"], clouds["iso_grav_anom"], footprint, known, transform
        )
        pixel_hash = canonical_pixel_sha256(prediction, ~footprint)
        stem = f"gemsdoe40-h4-euler-3d-dual-potential-kde-20261006-{pixel_hash[:8]}"
        output = args.candidate or (args.output_dir / f"{stem}.tif")
        write_candidate(
            output,
            prediction,
            args.sample,
            description=(
                "GEMSDOE40 H4 research-only candidate: depth-aware mutual-nearest SI=0 RTP magnetic "
                "and SI=-1 gravity Euler source-cloud KDE; proxy-gated, not organizer-scored"
            ),
        )
        receipt = validate_candidate(output, args.sample)
        if not receipt["valid"]:
            raise RuntimeError(f"candidate failed the exact sample-format contract: {receipt['issues']}")
        if receipt["canonical_pixels_sha256"] != pixel_hash:
            raise RuntimeError("canonical pixel hash changed during GeoTIFF round-trip")

    report = {
        "run_date_utc": "2026-10-06",
        "hypothesis": "H4 — depth-aware dual-potential-field Euler point-process consensus",
        "evidence_class": "candidate construction only; no proxy or holdout read during construction",
        "input_sha256": verified,
        "input_band_indexes": bands,
        "field_inputs": field_stats,
        "euler_windows": solve_stats,
        "method": h4_summary,
        "candidate_path": str(output),
        "candidate_sha256": receipt["sha256"],
        "canonical_pixels_sha256": pixel_hash,
        "format_receipt": receipt,
        "promotion_status": "NOT EVALUATED — do not submit until the separately locked holdout and uniqueness audit pass",
        "weekly_slot_used": False,
        "caveat": (
            "The SI=-1 gravity solve is an ordinary single-window Euler approximation to a finite gravity step; "
            "Reid & Thurston (2014) describe the fully generalized finite-step formulation as more complex. "
            "The candidate is research-only pending current-proxy holdout and refreshed prior-raster audit."
        ),
    }
    summary_path = args.summary or Path("docs/data/h4-build-20261006.json")
    write_json(summary_path, report)
    receipt_path = args.format_receipt or Path("docs/data/format-receipt-h4-20261006.json")
    write_json(receipt_path, receipt)
    return output, report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, default=Path("data/training_features.tif"))
    parser.add_argument("--sample", type=Path, default=Path("data/sample_submission.tif"))
    parser.add_argument("--labels", type=Path, default=Path("data/labels.tif"))
    parser.add_argument("--output-dir", type=Path, default=Path("docs/downloads"))
    parser.add_argument("--candidate", type=Path, default=None)
    parser.add_argument("--summary", type=Path, default=None)
    parser.add_argument("--format-receipt", type=Path, default=None)
    args = parser.parse_args()
    output, report = build(args)
    print(json.dumps({
        "candidate": str(output),
        "sha256": report["candidate_sha256"],
        "pixel_sha256": report["canonical_pixels_sha256"],
        "magnetic_solutions": report["method"]["magnetic"]["accepted_solutions"],
        "gravity_solutions": report["method"]["gravity"]["accepted_solutions"],
        "mutual_pairs": report["method"]["mutual_nearest_pairing"]["mutual_nearest_pairs"],
        "positive_pixels": report["method"]["kde"]["positive_in_footprint"],
        "format_valid": report["format_receipt"]["valid"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
