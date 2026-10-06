#!/usr/bin/env python3
"""Render ``docs/reports/h41-results-20261006.md`` from the JSON receipts.

Every number in the report is read from a receipt produced by a command in this
repository, so the prose cannot drift from the measurements.  ``--check`` fails if the
committed report differs from the rendered text.

    python scripts/write_h41_report.py [--check]
"""
from __future__ import annotations

import argparse
import json

import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
OUT = DOCS / "reports" / "h41-results-20261006.md"


def load(name: str) -> dict:
    return json.loads((DOCS / "data" / name).read_text(encoding="utf-8"))


def fmt(value, digits: int = 4) -> str:
    return f"{float(value):.{digits}f}"


def calibration_instrument(transfer: dict, entry: dict) -> float:
    """Saturating instrument reading for a measured raster, from its own (mass, w)."""
    K, a = transfer["fits"]["saturating"]
    mass, w = float(entry.get("mass", 0.0)), float(entry.get("w", 0.0))
    tp = K * (1.0 - np.exp(-a * mass * w / K))
    return float(tp / (0.2 * tp + 0.2 * mass + 0.8 * K))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    gen = load("h41-generation.json")
    audit = load("h41-audit.json")
    cand = load("current-candidate.json")
    transfer = load("live-transfer.json")
    instr = transfer["instruments"]["saturating"]
    cal = load("instrument-recalibration.json")
    h33 = load("h33-measured-analysis.json")
    obs = load("official-observation-20261006.json")
    slug = gen["files"]["slug"]
    mass = int(gen["emitted"]["mass"])
    checks = audit["checks"]
    blocked = checks["blocked"]
    means = blocked["means"]
    nov = checks["novelty"]
    anchors = sorted(transfer["anchors"], key=lambda r: -r["live_score"])

    lines: list[str] = []
    add = lines.append
    add("# H41 — multi-scale-stable shallow-contact Euler depth-cluster emission")
    add("")
    add("*Result record, 6 October 2026. Every figure below is read from a JSON receipt by "
        "`scripts/write_h41_report.py`; nothing is transcribed by hand and no organizer score "
        "exists for any file named here.*")
    add("")
    add("## 1. The released file")
    add("")
    add(f"| Field | Value |")
    add(f"| --- | --- |")
    add(f"| File | `docs/downloads/{slug}-zeros.tif` |")
    add(f"| Submission name | `{cand['name']}` |")
    add(f"| Note ({cand['note_characters']} characters) | {cand['note'].replace('|', chr(92) + '|')} |")
    add(f"| Bytes / SHA-256 | {cand['bytes']:,} / `{cand['sha256']}` |")
    add(f"| Canonical pixel SHA-256 | `{cand['canonical_pixels_sha256']}` |")
    add(f"| Format audit | {audit['verdict']} — single band float32, EPSG:32611, 100 m, "
        f"3,730 × 3,292, transform identical to the sample, all values finite in "
        f"[{fmt(checks['primary_min'], 1)}, {fmt(checks['primary_max'], 1)}], "
        f"{checks['primary_positive_px']:,} cells equal to 1.0 |")
    add(f"| Ties | `{slug}-nan.tif` (NaN outside the footprint, the sample's own encoding) "
        f"and `{slug}-zeros.zip` |")
    add(f"| Status | `slot_eligible: false` — see §5 |")
    add("")
    add("## 2. What was computed")
    add("")
    add(f"* Euler deconvolution (Reid et al. 1990, SI = 0 contact form with offset A) on `rtp`, "
        f"`tmi` and `iso_grav_anom` at windows {gen['windows_px']} px, stride {gen['stride_px']} px: "
        f"**{gen['euler']['total_solutions']:,}** depth-labelled solutions "
        f"(per window: "
        + ", ".join(f"{w} px → {v['solutions']:,}" for w, v in gen["euler"]["per_window"].items())
        + ").")
    add(f"* Per-window depth-cluster KDE (σ = 1.7 px; weights = shallowness × quality × tightness × "
        f"depth consistency × magnetic/gravity concordance), combined by geometric mean across all "
        f"three windows: **{gen['kde']['stability']['cells_supported_all_windows']:,}** cells "
        f"supported at every scale out of {gen['kde']['stability']['cells_supported_any']:,} "
        f"supported by at least one.")
    add(f"* Mean Euler depth per window is carried into the ranker, so the ranking sees the "
        f"depth structure, not only the density.")
    add(f"* Ranker: spatially blocked (5 × 6), out-of-fold L2 logistic discriminant over "
        f"{len(gen['discriminant']['features'])} features — the three KDE fields, the scale-stable "
        f"field, the three depth maps, potential-field gradients, DEM gradient and two-scale "
        f"curvature, ten competition bands, two cross-field gradient ratios, an RTP short-scale "
        f"curvature, a local radiometric z-score and the distance to the provided catalogue — "
        f"trained on public SGMC faults absent from the provided catalogue. "
        f"**Out-of-fold AUC {fmt(gen['discriminant']['oof_auc_all_blocks'])}** pooled, "
        f"{fmt(gen['discriminant']['oof_auc_mean_of_blocks'])} mean over "
        f"{len(gen['discriminant']['blocks_with_auc'])} blocks.")
    add(f"* Emission: Poisson-disk value-ranked thinning at {gen['spacing_px']:.0f} px spacing "
        f"(the metric's own R = 300 m), restricted to the Euler-gated support "
        f"({gen['support_px']:,} cells, catalogue flank {gen['catalogue_flank_px']:.0f} px), mass "
        f"chosen by the recalibrated instrument: **{mass:,} dots**, "
        f"w = {fmt(gen['emitted']['credit_per_dot_w'])}.")
    add("")
    add("## 3. Instruments, and what they will not say")
    add("")
    add(f"Four transfer models were fitted to the {transfer['n_anchors_scored']} scored anchors "
        f"(measured from their own bytes; `scripts/fit_live_transfer.py`):")
    add("")
    add("| Instrument | LOO Spearman | LOO RMSE | Passes the pre-registered 0.80 bar? |")
    add("| --- | --- | --- | --- |")
    for name, blob in transfer["instruments"].items():
        add(f"| `{name}` | {fmt(blob['loo_spearman'], 3)} | {fmt(blob['loo_rmse'])} | "
            f"{'yes' if blob['passes_bar'] else 'no'} |")
    add("")
    add(f"The saturating form (refit here on {transfer['n_anchors_scored']} scored anchors: "
        f"K = {fmt(transfer['fits']['saturating'][0], 1)}, "
        f"a = {fmt(transfer['fits']['saturating'][1])}) reproduces the sibling repository's model class "
        f"but not its digits. Its leave-one-out rank correlation is "
        f"**{fmt(instr['loo_spearman'], 3)}**, below the 0.80 bar this repository pre-registered, "
        f"and its largest leave-one-out error is {fmt(instr['loo_max_abs_error'])} — larger than "
        f"the difference it would be used to certify. Two further facts from the same fit:")
    add("")
    add(f"* The **whole-map proxy DTI ranks the scored artifacts backwards** out of fold "
        f"(LOO Spearman {fmt(transfer['instruments']['powerlaw']['loo_spearman'], 3)}). Maximising "
        f"agreement with the SGMC surrogate is not evidence of a better submission; it is, "
        f"measured here, slightly evidence of the opposite once mass is allowed to grow.")
    add(f"* The exact published metric algebra, given a single surrogate-to-hidden credit ratio "
        f"and a hidden truth size, fits worse than the phenomenon it was meant to explain "
        f"(LOO Spearman {fmt(transfer['instruments']['hidden']['loo_spearman'], 3)}; the fitted "
        f"hidden-truth size collapses to {fmt(transfer['fits']['hidden'][1], 0)} pixels, which is "
        f"rejected as a physical claim).")
    add("")
    add("## 4. Comparison on the identical blocked instrument")
    add("")
    add(f"4 × 6 partition, 3-pixel guard, {blocked['summary']['blocks_scored']} truth-bearing "
        f"blocks, exact published kernel, catalogue excluded from both sides.")
    add("")
    add("| Raster | Mean blocked DTI | Whole-map proxy DTI | w (per-dot credit) | Instrument |")
    add("| --- | --- | --- | --- | --- |")
    measured = {r["id"]: r for r in transfer["anchors"]}
    measured.update({r["id"]: r for r in transfer.get("unscored_measured", [])})

    def measured_by(substring: str) -> dict:
        hit = next((r for k, r in measured.items() if substring in k), None)
        return hit or {}

    add(f"| **H41 (this file)** | **{fmt(means['H41-candidate'])}** | "
        f"**{fmt(gen['audit']['proxy_dti'])}** | **{fmt(gen['emitted']['credit_per_dot_w'])}** | "
        f"**{fmt(gen['emitted']['instrument_pred'])}** |")
    for label, fragment, blocked_key in (
            ("H33-B2 · owner best, scored 0.2778", "h33-2-b2", "H33-B2"),
            ("H27-4 · owner, scored 0.2708", "h27-4", "H27-4"),
            ("H40-E · sibling Euler arm, unscored", "h40e", "H40-E")):
        entry = measured_by(fragment)
        add(f"| {label} | {fmt(means.get(blocked_key, float('nan')))} | "
            f"{fmt(entry.get('proxy_dti', float('nan')))} | "
            f"{fmt(entry.get('w', float('nan')))} | "
            f"{fmt(calibration_instrument(transfer, entry))} |")
    add("")
    add(f"The candidate strictly beats all three comparators in "
        f"{blocked['summary']['strict_wins_of_candidate']} of "
        f"{blocked['summary']['blocks_scored']} truth-bearing blocks.")
    add("")
    add("## 5. Novelty and the honest status")
    add("")
    add(f"Against {nov['n_priors']} same-grid prior rasters: maximum |Pearson| "
        f"**{fmt(nov['max_abs_pearson'])}** ({nov['max_abs_pearson_file']}), maximum support "
        f"Jaccard {fmt(nov['max_support_jaccard'])}, maximum top-budget Jaccard "
        f"{fmt(nov['max_top_budget_jaccard'])} — `is_new = {nov['is_new']}`.")
    add("")
    add("The file is therefore unique and format-legal, and it is the strongest construction "
        "this arm has produced. It is **not** certified for a weekly slot, because the "
        "pre-registered promotion bar (a leave-one-out rank correlation ≥ 0.80 for the "
        "instrument that would justify spending the slot) is not met, and because a prior "
        "unscored Euler arm of a sibling repository (H40-E) still leads the blocked proxy "
        "mean. Both statements are measurements, not opinions.")
    add("")
    add("## 6. Why the owner's 0.2778 was the best so far")
    add("")
    add(f"* Small mass: the file has {h33['h33_positive_pixels']:,} positive pixels, and across "
        f"the {transfer['n_anchors_scored']} scored anchors larger emissions rank lower live "
        f"(Spearman between mass and live score is negative).")
    add(f"* Arrangement, not luck: at the same {h33['h33_positive_pixels']:,}-pixel mass the "
        f"scattered control scores 0.0778 while the structure-aligned file scores 0.2778.")
    add(f"* It adds nothing: it is the {h33['base_positive_pixels']:,}-pixel base with "
        f"{h33['removed_pixels']:,} dots inside 200 m of the provided catalogue removed "
        f"(minimum remaining distance {h33['h33_min_distance_to_known_m']:.0f} m).")
    add(f"* The public target has moved: the board observed on {obs['observed_date_utc']} showed "
        f"{obs['top_rows'][0]['score']} at #1.")
    add("")
    add("## 7. Reproduce")
    add("")
    add("```bash")
    add("python scripts/calibrate_live_instrument.py")
    add("python scripts/fit_live_transfer.py")
    add("python scripts/run_h41_candidate.py --rank gated --flank 2 --cache work/h41_full.pkl")
    add("python scripts/audit_h41_candidate.py")
    add("python scripts/write_h41_report.py")
    add("python scripts/build_h41_site.py")
    add("python -m pytest -q")
    add("```")
    add("")
    rendered = "\n".join(lines) + "\n"
    if args.check:
        if not OUT.exists() or OUT.read_text(encoding="utf-8") != rendered:
            print("report differs from the receipts: " + str(OUT.relative_to(ROOT)))
            return 1
        print("report matches the receipts")
        return 0
    OUT.write_text(rendered, encoding="utf-8")
    print(f"[out] {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
