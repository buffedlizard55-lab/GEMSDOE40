#!/usr/bin/env python3
"""Validate H7 on the current SGMC proxy, all refreshed prior outputs, and the full pixel audit.

This script never uploads a file. It compares the candidate to every same-grid prior prediction
output (all prior output formats for holdout ranking), separately reports non-circular/format-valid
subsets, applies the frozen 24-block promotion rule against the highest prior score, and compares
raw/pixel patterns against all 306 hash-pinned prior blobs.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import rasterio

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gems40.pins import PINNED_FILES
from gemsdoe40.raster import validate_candidate, write_json
from gemsdoe40.research_holdout import (
    CURRENT_PROXY_SHA256,
    promotion_gate,
    read_proxy_truth,
    score_array_on_proxy,
)
from gemsdoe40.research_uniqueness import audit_candidate, file_sha256

TOP_BUDGET = 45_962
KNOWN_PROXY_CIRCULAR_BLOBS = {
    "9380203f3cb5edb9d096614f8dcf453e8f141b8f": "GEMSDOE3/7 SGMC-gap-only output, constructed from the same owner-derived SGMC family",
    "faa8ee8c0f7cb6b3a83385597c313e83a377c90c": "GEMSDOE4 proxy_only46 probability raster; proxy-derived research output",
}


def _is_proxy_circular(entry: dict[str, Any]) -> tuple[bool, str | None]:
    blob = entry["git_blob_sha"]
    if blob in KNOWN_PROXY_CIRCULAR_BLOBS:
        return True, KNOWN_PROXY_CIRCULAR_BLOBS[blob]
    classification = (entry.get("classification") or "").lower()
    if "circular" in classification:
        return True, classification
    for artifact in entry.get("artifacts", []):
        text = (artifact.get("path", "") + " " + artifact.get("repo", "")).lower()
        if "proxy_only" in text or "sgmc-gap" in text or "sgmc_gap" in text:
            return True, "filename/path explicitly identifies a proxy-only or SGMC-gap-derived output"
    return False, None


def _format_eligible(array: np.ndarray, dataset, footprint: np.ndarray, reference: dict[str, Any]) -> tuple[bool, list[str]]:
    issues: list[str] = []
    if dataset.dtypes[0] != "float32":
        issues.append("not float32")
    if dataset.nodata is None or not np.isnan(dataset.nodata):
        issues.append("nodata is not NaN")
    if not np.isfinite(array[footprint]).all():
        issues.append("non-finite in-footprint values")
    elif array[footprint].size and ((array[footprint] < 0).any() or (array[footprint] > 1).any()):
        issues.append("in-footprint value outside [0,1]")
    if np.isnan(array[~footprint]).sum() != int((~footprint).sum()):
        issues.append("outside-footprint cells are not all NaN")
    if (dataset.height, dataset.width) != (reference["height"], reference["width"]):
        issues.append("dimensions differ from template")
    if dataset.crs != reference["crs"] or dataset.transform != reference["transform"]:
        issues.append("CRS/transform differs from template")
    if dataset.count != 1:
        issues.append("not single band")
    return not issues, issues


def score_all_priors(
    inventory_path: Path,
    cache_dir: Path,
    truth: np.ndarray,
    valid: np.ndarray,
    labels: np.ndarray,
    reference: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    rows: list[dict[str, Any]] = []
    missing_cache: list[str] = []
    skipped_grid: list[dict[str, Any]] = []
    format_eligible_count = 0

    inventory_entries = inventory.get("unique_rasters", [])
    for index, entry in enumerate(inventory_entries, start=1):
        if index == 1 or index % 25 == 0 or index == len(inventory_entries):
            print(f"scoring prior TIFF {index}/{len(inventory_entries)}", flush=True)
        blob = entry["git_blob_sha"]
        cached = cache_dir / f"{blob}.tif"
        if not cached.exists():
            missing_cache.append(blob)
            continue
        raw_hash = file_sha256(cached)
        if entry.get("sha256") and raw_hash != entry["sha256"]:
            raise ValueError(f"cached prior SHA-256 mismatch for {blob}: {raw_hash} != {entry['sha256']}")
        with rasterio.open(cached) as dataset:
            same_grid = (
                dataset.count == 1
                and dataset.height == reference["height"]
                and dataset.width == reference["width"]
                and dataset.crs == reference["crs"]
                and dataset.transform == reference["transform"]
            )
            if not same_grid:
                skipped_grid.append({"git_blob_sha": blob, "sha256": raw_hash, "reason": "not a single-band exact-template grid"})
                continue
            array = dataset.read(1)
            eligible, format_issues = _format_eligible(array, dataset, valid, reference)
            score = score_array_on_proxy(array, truth, valid, labels)
        circular, circular_reason = _is_proxy_circular(entry)
        format_eligible_count += int(eligible)
        rows.append({
            "git_blob_sha": blob,
            "sha256": raw_hash,
            "bytes": cached.stat().st_size,
            "paths": entry.get("artifacts", []),
            "classification": entry.get("classification"),
            "proxy_circular": circular,
            "proxy_circular_reason": circular_reason,
            "same_grid": True,
            "format_eligible": eligible,
            "format_issues": format_issues,
            "in_footprint_nonzero": int(np.count_nonzero(array[valid])),
            "proxy_holdout": score,
        })

    scored = rows
    if not scored:
        raise RuntimeError("no same-grid prior output could be scored on the current proxy")
    scored.sort(key=lambda item: (item["proxy_holdout"]["pooled"]["score"], item["sha256"]), reverse=True)
    all_best = scored[0]
    noncircular = [row for row in scored if not row["proxy_circular"]]
    format_valid = [row for row in scored if row["format_eligible"]]
    combined_non_circular_format = [row for row in scored if row["format_eligible"] and not row["proxy_circular"]]
    summary = {
        "inventory_path": str(inventory_path),
        "inventory_unique_blobs": len(inventory.get("unique_rasters", [])),
        "inventory_current_head_repositories": inventory.get("current_head_refresh", {}).get("owner_public_repositories_matching_gemsdoe"),
        "inventory_current_head_tiff_paths": inventory.get("current_head_refresh", {}).get("current_head_tiff_paths_scanned"),
        "same_grid_prior_outputs_scored": len(rows),
        "prior_outputs_not_same_grid": len(skipped_grid),
        "prior_format_eligible_count": format_eligible_count,
        "proxy_circular_scored_count": sum(row["proxy_circular"] for row in scored),
        "missing_cache_count": len(missing_cache),
        "missing_cache_git_blobs": missing_cache,
        "highest_all_prior": _incumbent_summary(all_best),
        "highest_non_circular_prior": _incumbent_summary(noncircular[0]) if noncircular else None,
        "highest_format_eligible_prior": _incumbent_summary(format_valid[0]) if format_valid else None,
        "highest_non_circular_format_eligible_prior": _incumbent_summary(combined_non_circular_format[0]) if combined_non_circular_format else None,
        "top_10_all_priors": [_incumbent_summary(row) for row in scored[:10]],
        "top_10_non_circular": [_incumbent_summary(row) for row in noncircular[:10]],
        "not_same_grid": skipped_grid,
    }
    return rows, summary


def _incumbent_summary(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "git_blob_sha": row["git_blob_sha"],
        "sha256": row["sha256"],
        "score": row["proxy_holdout"]["pooled"]["score"],
        "positive_pixels": row["in_footprint_nonzero"],
        "format_eligible": row["format_eligible"],
        "format_issues": row["format_issues"],
        "proxy_circular": row["proxy_circular"],
        "proxy_circular_reason": row["proxy_circular_reason"],
        "paths": row["paths"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=Path, default=Path("docs/downloads/gemsdoe40-h7-rtp-euler-gravity-context-3d-kde-20261006-998f660f.tif"))
    parser.add_argument("--features", type=Path, default=Path("data/training_features.tif"))
    parser.add_argument("--sample", type=Path, default=Path("data/sample_submission.tif"))
    parser.add_argument("--labels", type=Path, default=Path("data/labels.tif"))
    parser.add_argument("--proxy", type=Path, default=Path("data/external/derived_sgmc_faults_100m_u8.tif"))
    parser.add_argument("--inventory", type=Path, default=Path("docs/data/prior_raster_inventory-20261006.json"))
    parser.add_argument("--cache", type=Path, default=Path("/tmp/gemsdoe40-prior-cache"))
    parser.add_argument("--report", type=Path, default=Path("docs/data/validation-h7-20261006.json"))
    parser.add_argument("--prior-scores-report", type=Path, default=Path("docs/data/prior-holdout-scores-h7-20261006.json"))
    parser.add_argument("--uniqueness-report", type=Path, default=Path("docs/data/uniqueness-audit-h7-20261006.json"))
    args = parser.parse_args()

    for key, path in (("features", args.features), ("sample", args.sample), ("labels", args.labels)):
        expected = PINNED_FILES[
            {"features": "training_features.tif", "sample": "example_submission.tif", "labels": "existing_faults.tif"}[key]
        ]["sha256"]
        actual = file_sha256(path)
        if actual != expected:
            raise ValueError(f"{key} SHA-256 mismatch: {actual} != {expected}")
    proxy_hash = file_sha256(args.proxy)
    if proxy_hash != CURRENT_PROXY_SHA256:
        raise ValueError(f"current proxy pin mismatch: {proxy_hash} != {CURRENT_PROXY_SHA256}")

    format_receipt = validate_candidate(args.candidate, args.sample)
    truth, valid, labels, proxy_info = read_proxy_truth(
        args.proxy,
        args.sample,
        args.labels,
        expected_proxy_sha256=CURRENT_PROXY_SHA256,
    )
    with rasterio.open(args.candidate) as candidate_ds:
        candidate = candidate_ds.read(1)
    candidate_score = score_array_on_proxy(candidate, truth, valid, labels)
    with rasterio.open(args.sample) as sample_ds:
        reference = {"height": sample_ds.height, "width": sample_ds.width, "crs": sample_ds.crs, "transform": sample_ds.transform}
    prior_rows, prior_summary = score_all_priors(
        args.inventory,
        args.cache,
        truth,
        valid,
        labels,
        reference,
    )
    write_json(args.prior_scores_report, {
        "proxy_sha256": proxy_hash,
        "proxy_status": proxy_info["proxy_status"],
        "candidate_sha256": format_receipt["sha256"],
        "inventory_path": str(args.inventory),
        "score_domain": "current pinned owner-derived proxy; 4x6 blocks with three-cell guards; exact labels==1 masking",
        "scores": prior_rows,
    })
    highest = max(prior_rows, key=lambda row: row["proxy_holdout"]["pooled"]["score"])
    gate = promotion_gate(candidate_score, highest["proxy_holdout"])
    uniqueness = audit_candidate(
        args.candidate,
        args.sample,
        args.inventory,
        args.cache,
        budget=TOP_BUDGET,
    )
    write_json(args.uniqueness_report, uniqueness)
    cache_complete = prior_summary["missing_cache_count"] == 0 and uniqueness["missing_prior_cache_count"] == 0
    slot_eligible = bool(gate["passed"] and format_receipt["valid"] and uniqueness["uniqueness_pass"] and cache_complete)
    decision = "ELIGIBLE FOR HUMAN REVIEW; NO SLOT USED" if slot_eligible else "HOLD — RESEARCH ONLY — DO NOT SUBMIT"
    report = {
        "run_date_utc": "2026-10-06",
        "hypothesis": "H7 — gravity-gradient-weighted magnetic Euler source cloud",
        "evidence_class": "owner-derived current SGMC proxy diagnostic; not organizer truth or an official score",
        "input_sha256": {
            "features": file_sha256(args.features),
            "sample": file_sha256(args.sample),
            "labels": file_sha256(args.labels),
            "current_proxy": proxy_hash,
            "candidate": format_receipt["sha256"],
        },
        "proxy_info": proxy_info,
        "candidate_path": str(args.candidate),
        "candidate_format_receipt": format_receipt,
        "candidate_proxy_holdout": candidate_score,
        "prior_corpus_holdout": prior_summary,
        "prior_holdout_scores_path": str(args.prior_scores_report),
        "promotion_gate_vs_highest_all_prior": gate,
        "uniqueness_audit_path": str(args.uniqueness_report),
        "uniqueness_summary": {
            "pass": uniqueness["uniqueness_pass"],
            "inventory_unique_blobs": uniqueness["inventory_unique_blobs"],
            "same_grid_comparisons": uniqueness["same_grid_comparisons"],
            "near_duplicate_count": uniqueness["near_duplicate_count"],
            "missing_prior_cache_count": uniqueness["missing_prior_cache_count"],
            "closest_prior_outputs": uniqueness["closest_prior_outputs"],
        },
        "cache_complete": cache_complete,
        "slot_eligible": slot_eligible,
        "final_action": decision,
        "weekly_slot_used": False,
        "organizer_score": None,
        "interpretation": "The all-prior incumbent is a local comparison only. Circular proxy-derived outputs are included in the conservative all-prior gate and shown separately; neither they nor the H7 score establish hidden/private-label performance.",
    }
    write_json(args.report, report)
    print(json.dumps({
        "candidate_score": candidate_score["pooled"]["score"],
        "candidate_blocks_won_vs_all_prior_best": gate["blocks_won"],
        "all_prior_incumbent": prior_summary["highest_all_prior"],
        "noncircular_incumbent": prior_summary["highest_non_circular_prior"],
        "same_grid_prior_outputs_scored": prior_summary["same_grid_prior_outputs_scored"],
        "format_eligible_prior_outputs": prior_summary["prior_format_eligible_count"],
        "uniqueness_pass": uniqueness["uniqueness_pass"],
        "closest_prior_outputs": uniqueness["closest_prior_outputs"][:3],
        "slot_eligible": slot_eligible,
        "final_action": decision,
        "report": str(args.report),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
