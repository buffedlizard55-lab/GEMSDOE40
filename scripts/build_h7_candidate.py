#!/usr/bin/env python3
"""Build the preregistered H7 magnetic Euler cloud with a bounded gravity context weight.

The H4 dual-Euler candidate was stopped before writing a TIFF because its locked raw-gravity SI=-1
solve had no depth-valid solutions. H7 is a separately registered fallback, not a relaxation of
H4: it uses only RTP for Euler source locations/depths and a continuous, non-zero gravity-gradient
percentile as a soft context multiplier. This construction stage does not open any proxy or prior
prediction raster.
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
    EulerCloud,
    potential_field_derivatives,
    project_weighted_cloud_to_kde,
    solve_euler_cloud,
    voxelize_and_weight,
)
from gemsdoe40.raster import band_index_by_name, canonical_pixel_sha256, validate_candidate, write_candidate, write_json
from gemsdoe40.research_euler import _gradient_components, fill_nearest


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
    verified = {}
    for key, path in paths.items():
        actual = file_sha256(path)
        expected = INPUT_SHA256[key]
        if actual != expected:
            raise ValueError(f"{key} SHA-256 mismatch: {actual} != pinned {expected}")
        verified[key] = actual
    return verified


def _align(features, sample, labels) -> None:
    ref = (sample.shape, sample.crs, sample.transform)
    if (features.shape, features.crs, features.transform) != ref:
        raise ValueError("feature stack is not exactly aligned with the sample grid")
    if (labels.shape, labels.crs, labels.transform) != ref:
        raise ValueError("known-label raster is not exactly aligned with the sample grid")
    if sample.crs is None or sample.crs.to_epsg() != 32611:
        raise ValueError(f"expected EPSG:32611 sample grid, got {sample.crs}")
    if abs(sample.transform.a - 100.0) > 1e-9 or abs(sample.transform.e + 100.0) > 1e-9:
        raise ValueError("sample grid must have 100 m north-up cells")


def gravity_context_factor(
    cloud,
    gradient_magnitude: np.ndarray,
    valid_footprint: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Evaluate a tie-averaged empirical percentile at each nearest source pixel.

    The percentile multiplier is locked to 0.5 + 0.5 * rank and therefore cannot zero a magnetic
    solution merely because its gravity context is weak.
    """
    h, w = valid_footprint.shape
    row = np.floor(cloud.row.astype(np.float64) + 0.5).astype(np.int64)
    col = np.floor(cloud.col.astype(np.float64) + 0.5).astype(np.int64)
    inside = (row >= 0) & (row < h) & (col >= 0) & (col < w)
    eligible = np.zeros(len(cloud), dtype=bool)
    eligible[inside] = valid_footprint[row[inside], col[inside]]
    population = np.asarray(gradient_magnitude, dtype=np.float32)[valid_footprint].astype(np.float64)
    if not np.isfinite(population).all() or population.size == 0:
        raise ValueError("gravity context gradient is not finite on its declared valid footprint")
    source_index = np.flatnonzero(eligible)
    if source_index.size == 0:
        raise RuntimeError("H7 stop: no magnetic source has an in-footprint gravity context pixel")
    values = np.asarray(gradient_magnitude, dtype=np.float32)[row[source_index], col[source_index]].astype(np.float64)
    if not np.isfinite(values).all():
        raise ValueError("gravity context gradient is not finite at retained Euler source pixels")
    ordered = np.sort(population)
    left = np.searchsorted(ordered, values, side="left")
    right = np.searchsorted(ordered, values, side="right")
    percentile_rank = (left.astype(np.float64) + right.astype(np.float64)) / (2.0 * ordered.size)
    factor = (0.5 + 0.5 * percentile_rank).astype(np.float32)
    return factor, eligible, {
        "gravity_gradient_magnitude_p05_p50_p95": np.percentile(population, [5, 50, 95]).tolist(),
        "source_gravity_gradient_percentile_p05_p50_p95": np.percentile(percentile_rank, [5, 50, 95]).tolist(),
        "source_context_factor_p05_p50_p95": np.percentile(factor, [5, 50, 95]).tolist(),
        "source_solutions_dropped_without_valid_gravity_context": int((~eligible).sum()),
        "source_solutions_retained_with_valid_gravity_context": int(eligible.sum()),
        "context_factor_rule": "0.5 + 0.5 * tie-averaged empirical CDF percentile; no hard threshold",
    }


def build(args: argparse.Namespace) -> tuple[Path, dict]:
    paths = {"features": args.features, "sample": args.sample, "labels": args.labels}
    verified = verify_inputs(paths)
    with rasterio.open(args.features) as features, rasterio.open(args.sample) as sample, rasterio.open(args.labels) as labels_ds:
        _align(features, sample, labels_ds)
        template = sample.read(1)
        footprint = np.isfinite(template)
        labels = labels_ds.read(1)
        known = (labels == 1) & footprint

        rtp_band = band_index_by_name(features, "rtp")
        gravity_band = band_index_by_name(features, "iso_grav_anom")
        rtp = features.read(rtp_band).astype(np.float32, copy=False)
        rtp_valid = footprint & np.isfinite(rtp) & (np.abs(rtp) < 1e30)
        filled_rtp, derivative_valid, mtx, mty, mtz = potential_field_derivatives(rtp)
        derivative_valid &= rtp_valid
        magnetic_cloud, solve_counts = solve_euler_cloud(
            filled_rtp,
            mtx,
            mty,
            mtz,
            derivative_valid,
            transform=features.transform,
            si=0.0,
            field="rtp",
        )
        if len(magnetic_cloud) == 0:
            raise RuntimeError(
                f"H7 stop: locked SI=0 magnetic Euler gates produced no accepted source points; {solve_counts}"
            )

        gravity_raw = features.read(gravity_band).astype(np.float32, copy=False)
        gravity_valid = footprint & np.isfinite(gravity_raw) & (np.abs(gravity_raw) < 1e30)
        gravity_filled, gravity_valid_from_data = fill_nearest(gravity_raw)
        gravity_valid &= gravity_valid_from_data
        gx, gy = _gradient_components(gravity_filled, resolution_m=100.0)
        gravity_gradient = np.hypot(gx, gy).astype(np.float32)
        context_factor, context_keep, context_summary = gravity_context_factor(
            magnetic_cloud,
            gravity_gradient,
            gravity_valid,
        )
        magnetic_cloud = EulerCloud(
            row=magnetic_cloud.row[context_keep],
            col=magnetic_cloud.col[context_keep],
            depth_m=magnetic_cloud.depth_m[context_keep],
            relative_residual=magnetic_cloud.relative_residual[context_keep],
            condition=magnetic_cloud.condition[context_keep],
            offset_px=magnetic_cloud.offset_px[context_keep],
            si=magnetic_cloud.si,
            field=magnetic_cloud.field,
        )
        voxels, voxel_stats = voxelize_and_weight(
            magnetic_cloud,
            footprint,
            features.transform,
            source_context_factor=context_factor,
        )
        support = voxels.weight > 0
        if not support.any():
            raise RuntimeError(
                "H7 stop: the locked 3-D magnetic source clusters have no positive weighted support; "
                "no threshold relaxation is permitted"
            )
        prediction, kde_stats = project_weighted_cloud_to_kde(
            voxels.row[support],
            voxels.col[support],
            voxels.depth_m[support],
            voxels.weight[support],
            footprint,
            known,
            stop_message="H7 stop: no positive gravity-context-weighted Euler support",
        )
        pixel_hash = canonical_pixel_sha256(prediction, ~footprint)
        stem = f"gemsdoe40-h7-rtp-euler-gravity-context-3d-kde-20261006-{pixel_hash[:8]}"
        output = args.candidate or args.output_dir / f"{stem}.tif"
        write_candidate(
            output,
            prediction,
            args.sample,
            description=(
                "GEMSDOE40 H7 research-only candidate: SI=0 RTP magnetic Euler source cloud, "
                "3-D depth clustering and bounded gravity-gradient context weight; not organizer-scored"
            ),
        )
        receipt = validate_candidate(output, args.sample)
        if not receipt["valid"]:
            raise RuntimeError(f"H7 failed exact sample-format validation: {receipt['issues']}")
        if receipt["canonical_pixels_sha256"] != pixel_hash:
            raise RuntimeError("canonical pixel hash changed during GeoTIFF round-trip")

    report = {
        "run_date_utc": "2026-10-06",
        "hypothesis": "H7 — gravity-gradient-weighted magnetic Euler source cloud",
        "evidence_class": "candidate construction only; no proxy or holdout read during construction",
        "input_sha256": verified,
        "layers": {"rtp": rtp_band, "iso_grav_anom": gravity_band},
        "euler_solve": {**magnetic_cloud.summary(), **solve_counts},
        "gravity_context": context_summary,
        "3d_voxel_cluster": voxel_stats,
        "kde": kde_stats,
        "candidate_path": str(output),
        "candidate_sha256": receipt["sha256"],
        "canonical_pixels_sha256": pixel_hash,
        "format_receipt": receipt,
        "promotion_status": "NOT EVALUATED — research-only until current-proxy holdout and refreshed uniqueness audit pass",
        "weekly_slot_used": False,
        "limitations": [
            "Gravity horizontal-gradient contrast is only a soft corroborator; lithologic contacts and basin-fill boundaries can create the same signature.",
            "The SI=-1 raw-gravity Euler arm H4 was infeasible under its preregistered condition and depth gates and is not represented in H7.",
            "No organizer score or hidden-test result is available for this locally built artifact.",
        ],
    }
    summary_path = args.summary or Path("docs/data/h7-build-20261006.json")
    write_json(summary_path, report)
    receipt_path = args.format_receipt or Path("docs/data/format-receipt-h7-20261006.json")
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
        "accepted_magnetic_solutions": report["euler_solve"]["accepted_solutions"],
        "clustered_voxels": report["3d_voxel_cluster"]["cluster_supported_voxels"],
        "positive_pixels": report["kde"]["positive_in_footprint"],
        "format_valid": report["format_receipt"]["valid"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
