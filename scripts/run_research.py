#!/usr/bin/env python3
"""Run a locked H1/H2/H2-B build, format check, uniqueness audit, and proxy holdout."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import rasterio

from gemsdoe40.euler import run_h1
from gemsdoe40.euler_h2 import run_h2
from gemsdoe40.holdout import promotion_gate, read_proxy_truth, score_array_on_proxy
from gemsdoe40.raster import validate_candidate, write_candidate, write_json
from gemsdoe40.uniqueness import audit_candidate, file_sha256


PINNED_INPUT_SHA256 = {
    "features": "4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5",
    "sample": "2176d08e485aa2cd2860ce8df539db4faf4d76163b38a4dd8c30a40454d35cbc",
    "labels": "7ba308ccdc4418b31a178f4f1ef21aaa6e152e4028f2f6f64b01f7eb25ae4093",
    "proxy": "26d142c4c93282cd94f6950ab96f22aeff59fbbea523d43d662e76fa1b161b5c",
}

KNOWN_PROXY_DERIVED_PRIORS = {
    "9380203f3cb5edb9d096614f8dcf453e8f141b8f": {
        "status": "known_proxy_circular",
        "basis": "The pinned artifact is named gapfinder-v2-sgmc-gap; its originating project identifies it as an SGMC-gap-only raster derived from USGS SGMC faults outside supplied labels.",
        "source_url": "https://github.com/buffedlizard55-lab/GEMSDOE3",
        "interpretation": "Retained in the preregistered all-prior incumbent rule, but not an independent generalization baseline on the owner-derived SGMC proxy.",
    }
}


def _metadata_ok(ds: rasterio.io.DatasetReader, reference: dict[str, Any]) -> bool:
    return (
        ds.count == 1
        and ds.height == reference["height"]
        and ds.width == reference["width"]
        and ds.crs == reference["crs"]
        and ds.transform == reference["transform"]
        and ds.dtypes[0] == "float32"
    )


def verify_pinned_inputs(paths: dict[str, Path]) -> dict[str, str]:
    """Fail before a long run if any external input differs from its recorded hash."""
    verified: dict[str, str] = {}
    for key, path in paths.items():
        digest = file_sha256(path)
        expected = PINNED_INPUT_SHA256[key]
        if digest != expected:
            raise ValueError(f"{key} input SHA-256 mismatch: {digest} != pinned {expected}")
        verified[key] = digest
    return verified


def _format_eligible(arr: np.ndarray, ds: rasterio.io.DatasetReader, footprint: np.ndarray) -> tuple[bool, list[str]]:
    issues: list[str] = []
    if arr.shape != footprint.shape:
        return False, ["shape mismatch"]
    if ds.nodata is None or not np.isnan(ds.nodata):
        issues.append("nodata is not NaN")
    if not np.isfinite(arr[footprint]).all():
        issues.append("non-finite prediction inside footprint")
    if arr[footprint].size and ((arr[footprint] < 0).any() or (arr[footprint] > 1).any()):
        issues.append("value outside [0,1] inside footprint")
    if np.isnan(arr[~footprint]).sum() != int((~footprint).sum()):
        issues.append("outside cells are not all NaN")
    return not issues, issues


def score_prior_inventory(
    inventory_path: Path,
    prior_cache: Path,
    truth: np.ndarray,
    valid: np.ndarray,
    labels: np.ndarray,
    reference: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    rows: list[dict[str, Any]] = []
    skipped_cache = []
    for entry in inventory.get("unique_rasters", []):
        blob = entry["git_blob_sha"]
        cached = prior_cache / f"{blob}.tif"
        if not cached.exists():
            skipped_cache.append(blob)
            continue
        with rasterio.open(cached) as ds:
            same_grid = _metadata_ok(ds, reference)
            if not same_grid:
                rows.append({
                    "git_blob_sha": blob,
                    "bytes": cached.stat().st_size,
                    "paths": entry.get("artifacts", []),
                    "same_grid": False,
                    "format_eligible": False,
                    "format_issues": ["not a single-band float32 exact-template grid"],
                    "scored": False,
                })
                continue
            arr = ds.read(1)
            format_ok, format_issues = _format_eligible(arr, ds, valid)
            raw_hash = file_sha256(cached)
            expected_hash = entry.get("sha256")
            if expected_hash and raw_hash != expected_hash:
                raise ValueError(f"prior raster SHA-256 mismatch for {blob}: {raw_hash} != {expected_hash}")
            scored = score_array_on_proxy(arr, truth, valid, labels)
            row = {
                "git_blob_sha": blob,
                "sha256": raw_hash,
                "bytes": cached.stat().st_size,
                "paths": entry.get("artifacts", []),
                "same_grid": True,
                "format_eligible": format_ok,
                "format_issues": format_issues,
                "in_footprint_nonzero": int(np.count_nonzero(arr[valid])),
                "scored": True,
                "proxy_holdout": scored,
            }
            rows.append(row)
    eligible = [row for row in rows if row.get("scored") and row.get("format_eligible")]
    if not eligible:
        raise RuntimeError("no format-eligible same-grid prior output could be scored")
    eligible.sort(key=lambda row: (row["proxy_holdout"]["pooled"]["score"], row["sha256"]), reverse=True)
    incumbent = eligible[0]
    summary = {
        "inventory_unique_rasters": len(inventory.get("unique_rasters", [])),
        "prior_cache_files_found": len(rows),
        "same_grid_scored": sum(bool(row.get("scored")) for row in rows),
        "format_eligible_and_scored": len(eligible),
        "skipped_missing_cache": len(skipped_cache),
        "skipped_cache_git_blob_sha": skipped_cache,
        "incumbent": {
            "git_blob_sha": incumbent["git_blob_sha"],
            "sha256": incumbent["sha256"],
            "score": incumbent["proxy_holdout"]["pooled"]["score"],
            "nonzero_pixels": incumbent["in_footprint_nonzero"],
            "paths": incumbent["paths"],
            "proxy_provenance": KNOWN_PROXY_DERIVED_PRIORS.get(
                incumbent["git_blob_sha"],
                {"status": "not_adjudicated", "interpretation": "Do not assume proxy independence without source review."},
            ),
        },
        "top_10_format_eligible": [
            {
                "git_blob_sha": row["git_blob_sha"],
                "sha256": row["sha256"],
                "score": row["proxy_holdout"]["pooled"]["score"],
                "nonzero_pixels": row["in_footprint_nonzero"],
                "paths": row["paths"],
                "proxy_provenance": KNOWN_PROXY_DERIVED_PRIORS.get(
                    row["git_blob_sha"],
                    {"status": "not_adjudicated", "interpretation": "Do not assume proxy independence without source review."},
                ),
            }
            for row in eligible[:10]
        ],
    }
    return rows, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", required=True, type=Path)
    parser.add_argument("--sample", required=True, type=Path)
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--proxy", required=True, type=Path)
    parser.add_argument("--prior-cache", required=True, type=Path)
    parser.add_argument("--hypothesis", choices=("H1", "H2", "H2B"), default="H1")
    parser.add_argument("--inventory", type=Path, default=Path("docs/data/prior_raster_inventory.json"))
    parser.add_argument("--candidate", type=Path, default=None)
    parser.add_argument("--validation-report", type=Path, default=Path("docs/data/validation.json"))
    parser.add_argument("--uniqueness-report", type=Path, default=Path("docs/data/uniqueness_audit.json"))
    parser.add_argument("--format-receipt", type=Path, default=Path("docs/data/format_receipt.json"))
    args = parser.parse_args()
    verified_input_hashes = verify_pinned_inputs({
        "features": args.features,
        "sample": args.sample,
        "labels": args.labels,
        "proxy": args.proxy,
    })

    if args.hypothesis == "H1":
        pred, method_summary = run_h1(args.features, args.sample, args.labels)
        default_candidate = Path("docs/downloads/gemsdoe40-h1-joint-euler-depth-kde-20261005.tif")
        description = "GEMSDOE40 preregistered H1 joint magnetic-gravity Euler depth-consensus KDE"
    elif args.hypothesis == "H2":
        pred, method_summary = run_h2(args.features, args.sample, args.labels)
        default_candidate = Path("docs/downloads/gemsdoe40-h2-tmi-euler-upward-persistence-20261005.tif")
        description = "GEMSDOE40 preregistered H2 TMI Euler upward-continuation persistence KDE"
    else:
        pred, method_summary = run_h2(args.features, args.sample, args.labels, support_only=True)
        default_candidate = Path("docs/downloads/gemsdoe40-h2b-tmi-euler-natural-support-20261005.tif")
        description = "GEMSDOE40 append-only H2-B TMI Euler upward-continuation persistence, all positive KDE support"
    args.candidate = args.candidate or default_candidate
    write_candidate(args.candidate, pred, args.sample, description=description)
    receipt = validate_candidate(args.candidate, args.sample)
    write_json(args.format_receipt, receipt)
    if not receipt["valid"]:
        raise RuntimeError(f"candidate failed exact sample-format validation: {receipt['issues']}")

    truth, valid, labels, proxy_info = read_proxy_truth(args.proxy, args.sample, args.labels)
    candidate_proxy = score_array_on_proxy(pred, truth, valid, labels)
    with rasterio.open(args.sample) as sample:
        reference = {"height": sample.height, "width": sample.width, "crs": sample.crs, "transform": sample.transform}
    prior_rows, prior_summary = score_prior_inventory(
        args.inventory, args.prior_cache, truth, valid, labels, reference
    )
    incumbent_row = next(row for row in prior_rows if row.get("sha256") == prior_summary["incumbent"]["sha256"])
    gate = promotion_gate(candidate_proxy, incumbent_row["proxy_holdout"])

    uniqueness = audit_candidate(args.candidate, args.sample, args.inventory, args.prior_cache)
    write_json(args.uniqueness_report, uniqueness)
    promotion = {
        **gate,
        "format_pass": bool(receipt["valid"]),
        "uniqueness_pass": bool(uniqueness["uniqueness_pass"]),
        "slot_eligible": bool(gate["passed"] and receipt["valid"] and uniqueness["uniqueness_pass"]),
        "final_action": "HOLD; DO NOT SUBMIT" if not (gate["passed"] and receipt["valid"] and uniqueness["uniqueness_pass"]) else "ELIGIBLE FOR USER REVIEW; no slot used by this process",
    }
    inputs = {
        "features_sha256": verified_input_hashes["features"],
        "sample_sha256": verified_input_hashes["sample"],
        "labels_sha256": verified_input_hashes["labels"],
        "sgmc_proxy_sha256": verified_input_hashes["proxy"],
        "prior_inventory_path": str(args.inventory),
    }
    report = {
        "run_date_utc": "2026-10-05",
        "evidence_class": "local SGMC-derived proxy holdout; not organizer truth or official score",
        "candidate_path": str(args.candidate),
        "candidate_sha256": receipt["sha256"],
        "candidate_canonical_pixels_sha256": receipt["canonical_pixels_sha256"],
        "inputs": inputs,
        "proxy": proxy_info,
        "proxy_independence_caveat": "Owner-derived SGMC mirror; not independently rebuilt or organizer truth. Some prior rasters are SGMC-derived. The frozen all-prior incumbent is retained as registered, but its score is not an independent generalization estimate.",
        "method_summary": method_summary,
        "prior_baselines": prior_summary,
        "candidate_proxy_holdout": candidate_proxy,
        "promotion_gate": promotion,
        "prior_raster_scores": prior_rows,
        "format_receipt_path": str(args.format_receipt),
        "uniqueness_report_path": str(args.uniqueness_report),
        "claims": {
            "official_leaderboard_score": None,
            "candidate_organizer_score": None,
            "weekly_submission_made": False,
            "preregistration_changed_after_result": False,
        },
    }
    write_json(args.validation_report, report)
    print(json.dumps({
        "candidate": str(args.candidate),
        "candidate_sha256": receipt["sha256"],
        "format_valid": receipt["valid"],
        "prior_rasters_scored": prior_summary["same_grid_scored"],
        "proxy_truth_pixels": proxy_info["off_catalogue_truth_pixels_after_exact_mask"],
        "incumbent": prior_summary["incumbent"],
        "candidate_proxy_dti": candidate_proxy["pooled"]["score"],
        "promotion_gate": promotion,
        "uniqueness_pass": uniqueness["uniqueness_pass"],
        "nearest_prior": uniqueness["closest_prior_outputs"][:3],
    }, indent=2))


if __name__ == "__main__":
    main()
