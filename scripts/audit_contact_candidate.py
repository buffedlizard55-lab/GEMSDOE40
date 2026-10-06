#!/usr/bin/env python3
"""Audit frozen H4 bytes; repeated scoring is reproduction, never a parameter search."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gemsdoe40.contact_audit import audit  # noqa: E402
from gemsdoe40.raster import validate_candidate, write_json  # noqa: E402
from gemsdoe40.research_holdout import read_proxy_truth, score_array_on_proxy, promotion_gate  # noqa: E402
from gemsdoe40.research_uniqueness import file_sha256  # noqa: E402
from acquire_data import PINS  # noqa: E402

BASELINES = {
    "frozen_sgmc_derived_incumbent": "9380203f3cb5edb9d096614f8dcf453e8f141b8f",
    "h33_b2_owner_reported_02778": "17a76895f68174cc93f3cb1597d25686f6d6bc67",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, default=ROOT / "work/h4")
    parser.add_argument("--inventory", type=Path, default=ROOT / "docs/data/prior-inventory-20261006.json")
    parser.add_argument("--prior-cache", type=Path, default=ROOT / "data/prior")
    args = parser.parse_args()
    args.work = args.work.resolve()
    generation = json.loads((args.work / "generation.json").read_text())
    if file_sha256(ROOT / "docs/research/h4-preregistration-20261006.md") != generation["preregistration_sha256"]:
        raise RuntimeError("frozen registration changed after generation")
    for name, expected in PINS.items():
        if file_sha256(ROOT / "data" / name) != expected:
            raise RuntimeError(f"pinned audit input changed: {name}")
    for source, expected in generation["code_sha256"].items():
        if file_sha256(ROOT / source) != expected:
            raise RuntimeError(f"generator code changed after the recorded run: {source}")
    candidate = ROOT / generation["candidate_format"]["path"]
    if file_sha256(candidate) != generation["candidate_format"]["sha256"]:
        raise RuntimeError("candidate bytes changed after generation")
    sample = ROOT / "data/sample_submission.tif"
    uniqueness = audit(candidate, sample, args.inventory, args.prior_cache)
    write_json(args.work / "uniqueness.json", uniqueness)
    if not uniqueness["uniqueness_pass"]:
        raise RuntimeError(f"refuse to publish/call new: incomplete corpus audit or near-duplicate; see {args.work / 'uniqueness.json'}")
    truth, valid, labels, proxy_info = read_proxy_truth(
        ROOT / "data/external/derived_sgmc_faults_100m_u8.tif", sample, ROOT / "data/labels.tif")
    with rasterio.open(candidate) as s:
        prediction = s.read(1)
    candidate_score = score_array_on_proxy(prediction, truth, valid, labels)
    baseline_scores = {}
    inventory = {e["git_blob_sha"]: e for e in json.loads(args.inventory.read_text())["unique_rasters"]}
    for name, blob in BASELINES.items():
        path = args.prior_cache / f"{blob}.tif"
        if file_sha256(path) != inventory[blob]["sha256"]:
            raise RuntimeError("baseline integrity mismatch")
        with rasterio.open(path) as s:
            prior = s.read(1)
        baseline_scores[name] = {"sha256": file_sha256(path), "artifact": inventory[blob]["artifacts"][0],
                                 "score": score_array_on_proxy(prior, truth, valid, labels)}
    prior_h2b = ROOT / "docs/downloads/gemsdoe40-h2b-tmi-euler-natural-support-20261005.tif"
    with rasterio.open(prior_h2b) as s:
        prior = s.read(1)
    baseline_scores["h2b_previous_research"] = {"sha256": file_sha256(prior_h2b),
                                               "score": score_array_on_proxy(prior, truth, valid, labels)}
    incumbent = baseline_scores["frozen_sgmc_derived_incumbent"]["score"]
    if not np.isclose(incumbent["pooled"]["score"], .8359066540883113, atol=1e-9, rtol=0):
        raise RuntimeError("frozen incumbent proxy replication failed")
    gate = promotion_gate(candidate_score, incumbent)
    nonempty = sum(b["n_truth"] > 0 for b in candidate_score["blocks"])
    gate["blocks_with_truth"] = nonempty
    gate["empty_blocks"] = len(candidate_score["blocks"]) - nonempty
    gate["strict_win_requirement_feasible"] = nonempty >= gate["blocks_required"]
    gate["irregularity"] = ("Historical frozen gate asks for 18/24 strict wins but only 16 blocks contain truth; "
                            "empty-block DTI is zero for every candidate. It cannot be passed. Preserved, not silently relaxed.")
    format_receipt = validate_candidate(candidate, sample)
    gate.update(format_pass=format_receipt["valid"], uniqueness_pass=uniqueness["uniqueness_pass"],
                slot_eligible=bool(gate["passed"] and gate["strict_win_requirement_feasible"] and format_receipt["valid"] and uniqueness["uniqueness_pass"]),
                weekly_submission_used=False, organizer_score=None)
    gate["final_action"] = "ELIGIBLE FOR REVIEW" if gate["slot_eligible"] else "HOLD — DO NOT SUBMIT"
    # Scores are now fixed. No ranking, window, threshold, or KDE tuning follows.
    validation = {"experiment": "H4", "evaluated_utc": datetime.now(timezone.utc).isoformat(),
                  "evidence_class": "Spatially blocked owner-derived SGMC proxy; not independent organizer truth",
                  "candidate_sha256": format_receipt["sha256"], "candidate_score": candidate_score,
                  "audit_code_sha256": {source: file_sha256(ROOT / source) for source in (
                      "scripts/audit_contact_candidate.py", "src/gemsdoe40/contact_audit.py",
                      "src/gemsdoe40/research_metric.py", "src/gemsdoe40/research_holdout.py", "src/gemsdoe40/raster.py")},
                  "preregistration_sha256": file_sha256(ROOT / "docs/research/h4-preregistration-20261006.md"),
                  "baselines": baseline_scores, "proxy": proxy_info, "promotion_gate": gate,
                  "circularity_warning": "Frozen incumbent used SGMC itself; its proxy score is circular, not a generalization estimate.",
                  "score_claim_policy": "0.2778 is owner-reported for H33-B2; public leaderboard rows identify accounts, not this TIFF."}
    if not format_receipt["valid"]:
        raise RuntimeError("final format reread failed")
    target = ROOT / "docs/downloads" / candidate.name
    if target.exists() and file_sha256(target) != file_sha256(candidate):
        raise RuntimeError("will not overwrite an existing named candidate with different bytes")
    cloud_source = args.work / "h4-euler-solutions.csv.gz"
    if file_sha256(cloud_source) != generation["cloud"]["sha256"]:
        raise RuntimeError("cloud bytes changed after generation")
    cloud_target = ROOT / "docs/downloads/h4-euler-solutions.csv.gz"
    if cloud_target.exists() and file_sha256(cloud_target) != file_sha256(cloud_source):
        raise RuntimeError("will not replace an existing named solution cloud with different bytes")
    # All byte, cloud, novelty and format checks precede publication.
    target.parent.mkdir(parents=True, exist_ok=True)
    for source, published in ((candidate, target), (cloud_source, cloud_target)):
        temporary = published.with_suffix(published.suffix + ".partial")
        shutil.copyfile(source, temporary)
        temporary.replace(published)
    format_receipt["path"] = str(target.relative_to(ROOT))
    name = "GEMSDOE40-H4-CONTACT-OFFSET-" + format_receipt["canonical_pixels_sha256"][:12]
    note = "H4 SI0+A rank-aware TMI+dGdz Euler; shallow cross-window depth-consensus KDE. New raw field; HOLD, proxy gate failed; no organizer score."
    assert len(note) <= 200
    release = {"experiment": "H4", "filename": target.name, "path": str(target.relative_to(ROOT)),
               "sha256": file_sha256(target), "bytes": target.stat().st_size,
               "canonical_pixels_sha256": format_receipt["canonical_pixels_sha256"],
               "name": name, "note": note, "note_characters": len(note),
               "status": gate["final_action"], "slot_eligible": gate["slot_eligible"],
               "organizer_score": None, "download_purpose": "Distinct submission-format research artifact; NOT approved for upload",
               "cloud": {**generation["cloud"], "path": str(cloud_target.relative_to(ROOT))}}
    generation["candidate_format"] = format_receipt
    generation["state"] = "FORMAT + RAW NOVELTY VERIFIED; " + gate["final_action"]
    for name_, obj in (("h4-generation.json", generation), ("h4-format.json", format_receipt),
                      ("h4-uniqueness.json", uniqueness), ("h4-validation.json", validation), ("current-candidate.json", release)):
        write_json(ROOT / "docs/data" / name_, obj)
    (ROOT / "docs/downloads/h4-submission-note.txt").write_text(name + "\n" + gate["final_action"] + "\n" + note + "\n")
    print(json.dumps({"release": release, "candidate_proxy_dti": candidate_score["pooled"]["score"],
                      "gate": gate}, indent=2), flush=True)


if __name__ == "__main__":
    main()
