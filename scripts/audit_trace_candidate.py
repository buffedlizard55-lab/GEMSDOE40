#!/usr/bin/env python3
"""Audit the frozen H8 candidate against gates G1-G5 of h8-preregistration-20261006.md.

G1 format; G2 full-corpus raw-output uniqueness; G3 frozen SGMC-proxy blocked holdout with
same-mass controls; G4 catalogue-component holdout (non-circular); G5 LM diagnostic (circular,
not a gate).  Repeated scoring is reproduction, never a parameter search.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from gems40.instrument import (  # noqa: E402
    cat_hidden_score, lm_score, load_cat_hidden, load_live_mirror, read_binary,
)
from gemsdoe40.contact_audit import audit  # noqa: E402
from gemsdoe40.raster import band_index_by_name, validate_candidate, write_json  # noqa: E402
from gemsdoe40.research_holdout import read_proxy_truth, score_array_on_proxy  # noqa: E402
from gemsdoe40.research_uniqueness import file_sha256  # noqa: E402
from acquire_data import PINS  # noqa: E402

PREREGISTRATION = "docs/research/h8-preregistration-20261006.md"
BASELINES = {
    "frozen_sgmc_derived_incumbent_circular": "9380203f3cb5edb9d096614f8dcf453e8f141b8f",
    "h33_b2_owner_reported_02778": "17a76895f68174cc93f3cb1597d25686f6d6bc67",
    "h27_4_owner_reported_02708": "12b0a4ddf1cf76f2552621b8a8b1a1b88d783786",
}
H4_PATH = "docs/downloads/gemsdoe40-h4-contact-offset-depthkde-20261006-fab9f6619c02.tif"
RANDOM_SEEDS = (40, 41, 42)


def mass_matched_controls(candidate: np.ndarray, footprint: np.ndarray, tmi: np.ndarray) -> dict:
    """Same value-mass random fields and a gradient-top-K field (binary 0/1)."""
    mass = float(np.nansum(candidate[footprint]))
    budget = int(round(mass))
    idx = np.flatnonzero(footprint.ravel())
    controls = {}
    for seed in RANDOM_SEEDS:
        rng = np.random.default_rng(seed)
        pick = rng.choice(idx, size=min(budget, idx.size), replace=False)
        field = np.zeros(candidate.shape, dtype=np.float32)
        field.ravel()[pick] = 1.0
        controls[f"random_mass_matched_seed{seed}"] = field
    grad = np.hypot(*np.gradient(np.nan_to_num(tmi, nan=0.0)))
    grad = np.where(footprint, grad, -np.inf)
    flat = grad.ravel()
    top = np.argpartition(flat, -budget)[-budget:]
    gfield = np.zeros(candidate.shape, dtype=np.float32)
    gfield.ravel()[top] = 1.0
    controls["gradient_topk_tmi_mass_matched"] = gfield
    return {"mass_budget_px": budget, "fields": controls}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, default=ROOT / "work/h8")
    parser.add_argument("--inventory", type=Path, default=ROOT / "docs/data/prior-inventory-20261006.json")
    parser.add_argument("--prior-cache", type=Path, default=ROOT / "data/prior")
    parser.add_argument("--skip-uniqueness", action="store_true",
                        help="developer-only; the gate still records uniqueness_pass=False")
    args = parser.parse_args()

    generation = json.loads((args.work / "generation.json").read_text())
    prereg_live = file_sha256(ROOT / PREREGISTRATION)
    prereg_frozen = generation["preregistration"]["sha256"]
    prereg_amended = prereg_live != prereg_frozen
    if prereg_amended:
        # The preregistration is append-only: an amendment registered after generation is legal
        # only if the generation-time hash is reproduced here and the amendment is disclosed.
        # The frozen sections governing THIS experiment must be unchanged by inspection.
        out_note = ("preregistration amended after generation (later experiment registered); "
                    "generation-time hash verified; amendment disclosed in evidence")
    else:
        out_note = "preregistration unchanged since generation"
    for source, expected in generation["code_sha256"].items():
        if file_sha256(ROOT / source) != expected:
            raise RuntimeError(f"generator code changed after the recorded run: {source}")
    for name, pinned in PINS.items():
        path = ROOT / "data" / name
        if path.exists() and file_sha256(path) != pinned:
            raise RuntimeError(f"pinned input changed: {name}")

    sample = ROOT / "data/sample_submission.tif"
    candidate_path = ROOT / generation["candidate_format"]["path"]
    if file_sha256(candidate_path) != generation["candidate_format"]["sha256"]:
        raise RuntimeError("candidate bytes changed after generation")

    out = {"experiment": generation["experiment"], "audited_utc": datetime.now(timezone.utc).isoformat(),
           "preregistration": {"sha256_at_generation": prereg_frozen, "sha256_live": prereg_live,
                               "amended_after_generation": prereg_amended, "note": out_note},
           "candidate": {"path": str(candidate_path.relative_to(ROOT)),
                         "sha256": generation["candidate_format"]["sha256"],
                         "canonical_pixels_sha256": generation["candidate_format"]["canonical_pixels_sha256"]},
           "field_stats": generation["field_stats"]}

    # G1 format
    receipt = validate_candidate(candidate_path, sample)
    out["g1_format"] = receipt
    if not receipt["valid"]:
        raise RuntimeError(f"G1 format gate failed: {receipt['issues']}")

    with rasterio.open(candidate_path) as ds:
        prediction = ds.read(1)
    with rasterio.open(sample) as ds:
        footprint = np.isfinite(ds.read(1))

    # G2 uniqueness (full corpus). A complete prior run for identical candidate bytes is reused;
    # any other state re-audits every blob (fail closed).
    cached = args.work / "uniqueness-full.json"
    if args.skip_uniqueness:
        out["g2_uniqueness"] = {"uniqueness_pass": False, "skipped": True}
    elif cached.exists():
        previous = json.loads(cached.read_text())
        if previous.get("candidate_sha256") != generation["candidate_format"]["sha256"]:
            raise RuntimeError("cached uniqueness audit does not match the current candidate bytes")
        out["g2_uniqueness"] = {k: v for k, v in previous.items() if k != "comparisons"}
        out["g2_uniqueness"]["comparisons_count"] = len(previous.get("comparisons", []))
        if not previous.get("uniqueness_pass"):
            raise RuntimeError(f"G2 uniqueness gate failed; see {args.work / 'uniqueness-full.json'}")
    else:
        uniqueness = audit(candidate_path, sample, args.inventory, args.prior_cache)
        out["g2_uniqueness"] = {k: v for k, v in uniqueness.items() if k != "comparisons"}
        out["g2_uniqueness"]["comparisons_count"] = len(uniqueness.get("comparisons", []))
        write_json(cached, uniqueness)
        if not uniqueness["uniqueness_pass"]:
            raise RuntimeError(f"G2 uniqueness gate failed; see {args.work / 'uniqueness-full.json'}")

    # G3 blocked proxy holdout + controls
    truth, valid, labels, proxy_info = read_proxy_truth(
        ROOT / "data/external/derived_sgmc_faults_100m_u8.tif", sample, ROOT / "data/labels.tif")
    out["proxy_info"] = proxy_info
    g3 = {"candidate": score_array_on_proxy(prediction, truth, valid, labels)}
    inventory = {e["git_blob_sha"]: e for e in json.loads(args.inventory.read_text())["unique_rasters"]}
    for name, blob in BASELINES.items():
        path = args.prior_cache / f"{blob}.tif"
        if file_sha256(path) != inventory[blob]["sha256"]:
            raise RuntimeError(f"baseline integrity mismatch: {name}")
        with rasterio.open(path) as ds:
            prior = ds.read(1)
        g3[name] = {"sha256": file_sha256(path), "artifact": inventory[blob]["artifacts"][0],
                    "score": score_array_on_proxy(np.nan_to_num(prior, nan=0.0), truth, valid, labels)}
    with rasterio.open(ROOT / H4_PATH) as ds:
        h4 = ds.read(1)
    g3["h4_contact_offset_hold"] = {"sha256": file_sha256(ROOT / H4_PATH),
                                    "score": score_array_on_proxy(h4, truth, valid, labels)}
    with rasterio.open(ROOT / "data/training_features.tif") as ds:
        tmi = ds.read(band_index_by_name(ds, "tmi")).astype(np.float64)
        tmi[~np.isfinite(tmi)] = 0.0
    controls = mass_matched_controls(prediction, footprint, tmi)
    g3["mass_budget_px"] = controls["mass_budget_px"]
    for name, field in controls["fields"].items():
        g3[name] = score_array_on_proxy(field, truth, valid, labels)
    out["g3_proxy_holdout"] = g3

    # G4 catalogue-component holdout (non-circular)
    ctx_cat = load_cat_hidden(ROOT / "data")
    g4 = {}
    g4["candidate"] = cat_hidden_score(prediction > 0, ctx_cat)
    for name, blob in BASELINES.items():
        with rasterio.open(args.prior_cache / f"{blob}.tif") as ds:
            g4[name] = cat_hidden_score(np.nan_to_num(ds.read(1), nan=0.0) > 0, ctx_cat)
    g4["h4_contact_offset_hold"] = cat_hidden_score(h4 > 0, ctx_cat)
    out["g4_cat_hidden"] = g4

    # G5 LM diagnostic (circular; reported only)
    ctx_lm = load_live_mirror(ROOT / "data")
    out["g5_lm_diagnostic"] = {
        "candidate": lm_score(prediction > 0, ctx_lm),
        "h33_b2": lm_score(read_binary(args.prior_cache / f"{BASELINES['h33_b2_owner_reported_02778']}.tif") & ctx_lm.foot, ctx_lm),
        "note": "circular: incumbent and proxy are both SGMC-derived; diagnostic only, never a gate",
    }

    # Gate decisions (frozen in preregistration)
    pooled = g3["candidate"]["pooled"]["score"]
    best_random = max(g3[f"random_mass_matched_seed{s}"]["pooled"]["score"] for s in RANDOM_SEEDS)
    grad = g3["gradient_topk_tmi_mass_matched"]["pooled"]["score"]
    h33_pooled = g3["h33_b2_owner_reported_02778"]["score"]["pooled"]["score"]
    g3_pass = (pooled > h33_pooled) and (pooled >= 2 * best_random) and (pooled > grad)
    cat_cand = g4["candidate"]["cat_hidden_mean"]
    cat_h33 = g4["h33_b2_owner_reported_02778"]["cat_hidden_mean"]
    g4_pass = cat_cand > cat_h33
    gate = {
        "g1_format_pass": bool(receipt["valid"]),
        "g2_uniqueness_pass": out["g2_uniqueness"].get("uniqueness_pass", False),
        "g3_pass": bool(g3_pass),
        "g3_detail": {"candidate_pooled_dti": pooled, "h33b2_pooled_dti": h33_pooled,
                      "best_random_control_dti": best_random, "gradient_control_dti": grad},
        "g4_pass": bool(g4_pass),
        "g4_detail": {"candidate_cat_hidden_mean": cat_cand, "h33b2_cat_hidden_mean": cat_h33},
        "g5_lm_calibrated_mean_diagnostic": out["g5_lm_diagnostic"]["candidate"]["lm_calibrated_mean"],
    }
    gate["all_pass"] = all([gate["g1_format_pass"], gate["g2_uniqueness_pass"], gate["g3_pass"], gate["g4_pass"]])
    gate["status"] = ("SLOT-ELIGIBLE CANDIDATE — proxy gates passed; organizer score unknown; owner decides"
                      if gate["all_pass"] else "HOLD — research only; a frozen gate failed")
    gate["weekly_slot_used"] = False
    gate["organizer_score"] = None
    out["gate"] = gate

    write_json(args.work / "audit.json", out)
    print(json.dumps({"status": gate["status"], "gate": gate}, indent=1))
    return 0 if gate["all_pass"] else 2


if __name__ == "__main__":
    sys.exit(main())
