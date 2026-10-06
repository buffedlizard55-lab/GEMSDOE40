#!/usr/bin/env python3
"""Validate the frozen H8 artifact: integrity, format, uniqueness, and the
registered replacement gate (proxy diagnostics + quadrant catalogue skill).

The verdict is computed, never negotiated: the gate rule was registered in
docs/research/h8-preregistration-20261006.md before any score was read.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from gemsdoe40.contact_audit import audit  # noqa: E402
from gemsdoe40.raster import validate_candidate, write_json  # noqa: E402
from gemsdoe40.research_holdout import read_proxy_truth, score_array_on_proxy  # noqa: E402
from gemsdoe40.research_uniqueness import file_sha256  # noqa: E402
from gemsdoe40 import holdout  # noqa: E402
from acquire_data import PINS  # noqa: E402

H33B2_BLOB = "17a76895f68174cc93f3cb1597d25686f6d6bc67"
INCUMBENT_BLOB = "9380203f3cb5edb9d096614f8dcf453e8f141b8f"
H33B2_POOLED = 0.091550266181  # measured on this same frozen instrument (H4 validation)
BLOCK_WIN_REQUIREMENT = 10


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, default=ROOT / "work/h8")
    parser.add_argument("--inventory", type=Path, default=ROOT / "docs/data/prior-inventory-20261006.json")
    parser.add_argument("--prior-cache", type=Path, default=ROOT / "data/prior")
    parser.add_argument("--skip-uniqueness", action="store_true")
    args = parser.parse_args()

    generation = json.loads((args.work / "generation.json").read_text())
    # integrity: preregistration and generator code unchanged since the run
    if file_sha256(ROOT / generation["config"]["preregistration"]) != generation["preregistration_sha256"]:
        raise RuntimeError("preregistration changed after generation")
    for source, expected in generation["code_sha256"].items():
        if file_sha256(ROOT / source) != expected:
            raise RuntimeError(f"generator code changed after the recorded run: {source}")
    for rel, expected in PINS.items():
        if (ROOT / "data" / rel).exists() and file_sha256(ROOT / "data" / rel) != expected:
            raise RuntimeError(f"pinned input changed: {rel}")

    nan_tif = ROOT / generation["output"]["nan_tif"]["path"]
    zeros_tif = ROOT / generation["output"]["zeros_tif"]["path"]
    sample = ROOT / "data/sample_submission.tif"

    result: dict = {"generation_basename": generation["basename"], "checks": {}}

    # 1. format receipts (independent re-open)
    result["checks"]["format_nan"] = validate_candidate(nan_tif, sample)
    with rasterio.open(zeros_tif) as z:
        arr = z.read(1)
        finite = np.isfinite(arr).all()
        in01 = bool((arr[np.isfinite(arr)].min() >= 0.0) and (arr.max() <= 1.0))
        nodata_unset = z.nodata is None
        with rasterio.open(sample) as s:
            footprint = np.isfinite(s.read(1))
            zeros_outside = bool((arr[~footprint] == 0).all())
            shape_ok = arr.shape == footprint.shape
            crs_ok = z.crs == s.crs
            tr_ok = z.transform == s.transform
        result["checks"]["format_zeros_portal"] = {
            "path": str(zeros_tif), "dtype": z.dtypes[0], "all_finite": bool(finite),
            "in_0_1": in01, "nodata_unset": bool(nodata_unset), "zeros_outside": zeros_outside,
            "shape_ok": shape_ok, "crs_ok": bool(crs_ok), "transform_ok": bool(tr_ok),
            "valid": bool(finite and in01 and nodata_unset and zeros_outside and shape_ok and crs_ok and tr_ok),
        }
    if not result["checks"]["format_nan"]["valid"] or not result["checks"]["format_zeros_portal"]["valid"]:
        raise RuntimeError("format receipt failed; refusing to score")

    result["candidate_sha256"] = file_sha256(nan_tif)
    result["candidate_canonical_sha256"] = result["checks"]["format_nan"]["canonical_pixels_sha256"]
    result["candidate_format"] = {
        "valid": bool(result["checks"]["format_nan"]["valid"]),
        "sha256": result["candidate_sha256"],
        "bytes": nan_tif.stat().st_size,
    }

    # 2. synthetic suite is pytest-based; re-run is cheap and recorded by CI.
    result["checks"]["synthetic_suite"] = "tests/test_asa_spi.py (S1-S5), run in CI; see pytest receipt"

    # 3. uniqueness against the dated frozen corpus
    if args.skip_uniqueness:
        result["checks"]["uniqueness"] = {"skipped": True}
    else:
        uniqueness = audit(nan_tif, sample, args.inventory, args.prior_cache)
        write_json(args.work / "uniqueness.json", uniqueness)
        result["checks"]["uniqueness"] = {
            "pass": uniqueness["uniqueness_pass"],
            "inventory_blobs": uniqueness.get("inventory_blobs"),
            "hashed_and_audited_blobs": uniqueness.get("hashed_and_audited_blobs"),
            "full_size_raw_comparisons": uniqueness.get("full_size_raw_comparisons"),
            "max_abs_pearson": uniqueness.get("max_abs_pearson"),
            "max_top_jaccard": uniqueness.get("max_top_jaccard"),
            "max_top_containment": uniqueness.get("max_top_containment"),
            "detail_path": "work/h8/uniqueness.json",
        }
        if not uniqueness["uniqueness_pass"]:
            result["verdict"] = "REJECT_AS_NEW: near-duplicate or incomplete corpus audit"
            write_json(args.work / "validation.json", result)
            print(json.dumps(result["verdict"]))
            return

    # 4. proxy diagnostics + registered replacement gate
    truth, valid, labels, proxy_info = read_proxy_truth(
        ROOT / "data/external/derived_sgmc_faults_100m_u8.tif", sample, ROOT / "data/labels.tif")
    with rasterio.open(nan_tif) as s:
        prediction = np.nan_to_num(s.read(1), nan=0.0)
    candidate_score = score_array_on_proxy(prediction, truth, valid, labels)

    inventory = {e["git_blob_sha"]: e for e in json.loads(args.inventory.read_text())["unique_rasters"]}
    comparators = {}
    with rasterio.open(args.prior_cache / f"{H33B2_BLOB}.tif") as s:
        h33 = np.nan_to_num(s.read(1), nan=0.0)
    comparators["h33_b2"] = score_array_on_proxy(h33, truth, valid, labels)
    with rasterio.open(args.prior_cache / f"{INCUMBENT_BLOB}.tif") as s:
        incumbent = np.nan_to_num(s.read(1), nan=0.0)
    comparators["frozen_sgmc_incumbent_circular"] = score_array_on_proxy(incumbent, truth, valid, labels)
    h4_path = ROOT / "docs/downloads/gemsdoe40-h4-contact-offset-depthkde-20261006-fab9f6619c02.tif"
    if h4_path.exists():
        with rasterio.open(h4_path) as s:
            h4 = np.nan_to_num(s.read(1), nan=0.0)
        comparators["h4_previous_session"] = score_array_on_proxy(h4, truth, valid, labels)

    cand_pooled = float(candidate_score["pooled"]["score"])
    h33_pooled = float(comparators["h33_b2"]["pooled"]["score"])
    if not np.isclose(h33_pooled, H33B2_POOLED, atol=1e-6):
        raise RuntimeError(f"H33-B2 proxy replication failed: {h33_pooled} != {H33B2_POOLED}")
    truth_blocks = [(i, b, cb) for i, (b, cb) in enumerate(
        zip(comparators["h33_b2"]["blocks"], candidate_score["blocks"])) if b["n_truth"] > 0]
    wins = sum(1 for _, b, cb in truth_blocks if cb["score"] > b["score"])
    pooled_ok = cand_pooled > H33B2_POOLED
    wins_ok = wins >= BLOCK_WIN_REQUIREMENT
    gate = {
        "rule": "registered replacement gate (h8-preregistration-20261006.md)",
        "candidate_pooled_proxy_dti": cand_pooled,
        "h33_b2_pooled_proxy_dti": h33_pooled,
        "pooled_threshold": H33B2_POOLED,
        "pooled_pass": bool(pooled_ok),
        "truth_bearing_blocks": len(truth_blocks),
        "block_wins_vs_h33_b2": int(wins),
        "block_win_requirement": BLOCK_WIN_REQUIREMENT,
        "block_pass": bool(wins_ok),
        "gate_pass": bool(pooled_ok and wins_ok),
        "incumbent_pooled_proxy_dti_circular": float(
            comparators["frozen_sgmc_incumbent_circular"]["pooled"]["score"]),
    }

    # 5. quadrant catalogue skill vs controls (diagnostic, off-catalogue target caveat applies)
    with rasterio.open(sample) as s:
        fp = np.isfinite(s.read(1))
    labels_bool = (labels == 1) & fp
    quad = holdout.blocked_dti(prediction, labels_bool, fp)
    rnd = holdout.same_mass_random((prediction > 0).astype(np.float32), fp, np.random.default_rng(40))
    quad_random = holdout.blocked_dti(rnd, labels_bool, fp)

    result["proxy_diagnostics"] = {
        "candidate_pooled": candidate_score["pooled"],
        "comparators": {k: v["pooled"] for k, v in comparators.items()},
        "proxy_info": proxy_info,
    }
    result["gate"] = gate
    result["catalogue_skill_diagnostic"] = {
        "candidate_mean_dti": quad["mean_dti"],
        "random_same_mass_mean_dti": quad_random["mean_dti"],
        "note": "catalogue is the KNOWN mask, not the hidden truth; this measures fault-finding style only",
    }
    result["verdict"] = ("ELIGIBLE_FOR_SUBMISSION_REVIEW (human portal decision required)"
                         if gate["gate_pass"] else "HOLD — DO NOT SUBMIT (registered gate failed)")
    result["promotion_gate"] = {
        "rule": gate["rule"],
        "passed": bool(gate["gate_pass"]),
        "slot_eligible": bool(gate["gate_pass"]),
        "organizer_score": None,
        "weekly_submission_used": False,
        "final_action": "ELIGIBLE_FOR_SUBMISSION_REVIEW" if gate["gate_pass"] else "HOLD — DO NOT SUBMIT",
        "candidate_pooled_proxy_dti": gate["candidate_pooled_proxy_dti"],
        "pooled_threshold": gate["pooled_threshold"],
        "block_wins_vs_h33_b2": gate["block_wins_vs_h33_b2"],
        "block_win_requirement": gate["block_win_requirement"],
    }
    write_json(args.work / "validation.json", result)
    write_json(ROOT / "docs/data/h8-validation.json", result)
    print(json.dumps({"verdict": result["verdict"], "gate": gate}, indent=2))


if __name__ == "__main__":
    main()
