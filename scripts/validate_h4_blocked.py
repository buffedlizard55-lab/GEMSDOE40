#!/usr/bin/env python3
"""Validate the H4 candidates on the repository's frozen blocked instrument.

Stage 1 re-verifies the instrument itself against the published (owner-reported)
score orderings; if the instrument cannot reproduce them, no decision may be
taken from it.  Stage 2 scores the shipped H4 raster and its reserve variants on
the four spatially blocked cells of the LM instrument and on the conventional
catalogue-component holdout, and reports per-fold wins against the best prior.

Writes ``docs/data/h4_blocked_validation.json``.
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

from gems40.instrument import cat_hidden_score, lm_score, load_cat_hidden, load_live_mirror, read_binary  # noqa: E402
from gemsdoe40.corpus import PRIORS, local_path  # noqa: E402

# (label, local corpus key, owner-reported score) - the published live chain used
# by the sibling repository to calibrate this instrument.
LIVE_CHAIN = (
    ("h33-2-b2", "h33-2-b2", 0.2778),
    ("h27-4-r1-solo-d2-8", "h27-4-r1-solo-d2-8", 0.2708),
    ("h32-1-prethin-tip-euler-d2-8", "h32-1-prethin-tip-euler-d2-8", 0.2649),
    ("dotted-h19-5-d2-8", "dotted-h19-5-d2-8", 0.2600),
    ("dotted-h19-5-d1-5", "dotted-h19-5-d1-5", 0.2477),
    ("topo-gap-closure-v2-d1-5", "topo-gap-closure-v2-d1-5", 0.2449),
    ("h19-5", "h19-5", 0.1922),
    ("h16-1-topo-geophys-ridges", "h16-1-topo-geophys-ridges", 0.1855),
    ("h28-dotted-ridge", "h28-dotted-ridge", 0.1839),
)
ORDERINGS = (
    ("h33-2-b2", "h27-4-r1-solo-d2-8"), ("h27-4-r1-solo-d2-8", "dotted-h19-5-d2-8"),
    ("dotted-h19-5-d2-8", "dotted-h19-5-d1-5"), ("dotted-h19-5-d1-5", "h19-5"),
    ("h19-5", "h16-1-topo-geophys-ridges"), ("h16-1-topo-geophys-ridges", "h28-dotted-ridge"),
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=str(ROOT / "data"))
    ap.add_argument("--candidates", nargs="*", default=[],
                    help="candidate GeoTIFFs; defaults to the shipped H4 file in docs/downloads")
    ap.add_argument("--out", default=str(ROOT / "docs" / "data" / "h4_blocked_validation.json"))
    args = ap.parse_args()

    ctx = load_live_mirror(args.data)
    cat = load_cat_hidden(args.data)
    print(f"[LM] 4 blocked cells · off-catalogue SGMC truth px = "
          f"{sum(int(c.truth.sum()) for c in ctx.cells)}", flush=True)

    report = dict(instrument=dict(
        kind="LM (spatially blocked quadrants, off-catalogue SGMC truth, prevalence calibrated)",
        cells={c.key: dict(domain_px=int(c.domain.sum()), truth_px=int(c.truth.sum()))
               for c in ctx.cells},
        catalogue_hidden_cells={c.key: dict(truth_px=int(c.truth.sum())) for c in cat.cells},
    ), stage1={}, candidates={})

    scores: dict[str, dict] = {}
    for label, key, lb in LIVE_CHAIN:
        path = local_path(next(p for p in PRIORS if p.key == key))
        if not path.exists():
            print(f"  ! missing {path}", flush=True)
            continue
        m = read_binary(path)
        r = lm_score(m, ctx)
        r["leaderboard_dti"] = lb
        scores[label] = r
    print(f"\n{'reference (owner-reported live)':34s} {'live':>7s} {'LM':>8s} {'LM-cal':>8s} {'N_px':>7s}",
          flush=True)
    for label, _, lb in LIVE_CHAIN:
        if label in scores:
            r = scores[label]
            print(f"{label:34s} {lb:7.4f} {r['lm_mean']:8.5f} {r['lm_calibrated_mean']:8.5f} "
                  f"{r['emitted_pixels']:7d}", flush=True)

    checks = []
    for hi, lo in ORDERINGS:
        if hi not in scores or lo not in scores:
            continue
        checks.append(dict(
            ordering=f"{hi} > {lo}",
            reproduced_raw=bool(scores[hi]["lm_mean"] > scores[lo]["lm_mean"]),
            reproduced_calibrated=bool(scores[hi]["lm_calibrated_mean"] >
                                       scores[lo]["lm_calibrated_mean"]),
            lm_margin=float(scores[hi]["lm_mean"] - scores[lo]["lm_mean"]),
            lm_calibrated_margin=float(scores[hi]["lm_calibrated_mean"] -
                                       scores[lo]["lm_calibrated_mean"])))
    n_ok = sum(c["reproduced_calibrated"] for c in checks)
    n_ok_raw = sum(c["reproduced_raw"] for c in checks)
    report["stage1"] = dict(orderings=checks, reproduced_calibrated=f"{n_ok}/{len(checks)}",
                            reproduced_raw=f"{n_ok_raw}/{len(checks)}",
                            reference_scores={k: dict(lm=v["lm_mean"],
                                                      lm_calibrated=v["lm_calibrated_mean"],
                                                      live=v["leaderboard_dti"])
                                              for k, v in scores.items()})
    print(f"\n[stage 1] orderings reproduced — LM {n_ok_raw}/{len(checks)}, "
          f"LM-calibrated {n_ok}/{len(checks)}", flush=True)

    # ---- calibration of the instrument onto the public leaderboard scale ----
    pts = [(v["lm_calibrated_mean"], v["leaderboard_dti"]) for v in scores.values()]
    xs = np.array([p[0] for p in pts]); ys = np.array([p[1] for p in pts])
    A = np.vstack([xs, np.ones_like(xs)]).T
    slope, intercept = np.linalg.lstsq(A, ys, rcond=None)[0]
    pred = slope * xs + intercept
    ss_res = float(((ys - pred) ** 2).sum()); ss_tot = float(((ys - ys.mean()) ** 2).sum())
    ratios = ys / xs
    report["stage1"]["leaderboard_projection"] = dict(
        n=len(pts), slope=float(slope), intercept=float(intercept),
        r2=1.0 - ss_res / ss_tot, rmse=float(np.sqrt(ss_res / len(ys))),
        mean_ratio=float(ratios.mean()), sd_ratio=float(ratios.std(ddof=1)),
        note=("LM-calibrated score -> owner-reported leaderboard DTI. All owner-reported "
              "scores are assumptions (the public leaderboard exposes no filenames/hashes)."))
    print(f"\n[stage 2] projection to leaderboard scale: DTI = {slope:.4f}*LMcal "
          f"{intercept:+.4f}  R2={1.0 - ss_res / ss_tot:.3f}  rmse={np.sqrt(ss_res / len(ys)):.4f}  "
          f"mean ratio={ratios.mean():.3f}", flush=True)
    print(f"          a 0.3262 leaderboard entry implies LM-cal ~= {(0.3262 - intercept) / slope:.4f}",
          flush=True)

    # reference: the best prior by the calibrated instrument
    ref_key, ref = max(scores.items(), key=lambda kv: kv[1]["lm_calibrated_mean"])
    candidates = args.candidates or [str(p) for p in sorted((ROOT / "docs" / "downloads").glob(
        "*euler-line-ring-pruned*-zeros.tif"))]
    for path in candidates:
        p = Path(path)
        if not p.exists():
            print(f"  ! missing candidate {p}", flush=True)
            continue
        m = read_binary(p)
        r = lm_score(m, ctx)
        r["cat_hidden"] = cat_hidden_score(m, cat)["cat_hidden_mean"]
        wins = sum(1 for k in r["lm_calibrated_per_fold"]
                   if r["lm_calibrated_per_fold"][k] >
                   scores[ref_key]["lm_calibrated_per_fold"][k])
        r["projected_leaderboard_dti"] = float(slope * r["lm_calibrated_mean"] + intercept)
        r["projected_leaderboard_dti_ci"] = [float(max(0.0, (slope * r["lm_calibrated_mean"] + intercept)
                                                       - 2 * np.sqrt(ss_res / len(ys)))),
                                             float(slope * r["lm_calibrated_mean"] + intercept
                                                   + 2 * np.sqrt(ss_res / len(ys)))]
        r["reference"] = ref_key
        r["folds_better_than_reference"] = wins
        r["beats_reference"] = bool(wins >= 3)
        report["candidates"][p.name] = r
        print(f"\n[{p.name}]", flush=True)
        print(f"  LM {r['lm_mean']:.5f}  LM-cal {r['lm_calibrated_mean']:.5f}  "
              f"CAT-HID {r['cat_hidden']:.5f}  N_px {r['emitted_pixels']}  "
              f"on_catalogue {r['on_catalogue_pixels']}", flush=True)
        print(f"  reference {ref_key}: LM-cal {scores[ref_key]['lm_calibrated_mean']:.5f}; "
              f"folds better: {wins}/4", flush=True)
        for k in sorted(r["lm_calibrated_per_fold"]):
            print(f"    {k}: candidate {r['lm_calibrated_per_fold'][k]:.5f} vs "
                  f"reference {scores[ref_key]['lm_calibrated_per_fold'][k]:.5f}", flush=True)
    Path(args.out).write_text(json.dumps(report, indent=1))
    print(f"\nwrote {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
