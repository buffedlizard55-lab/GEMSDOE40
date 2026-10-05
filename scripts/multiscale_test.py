#!/usr/bin/env python
"""H40-4: does a multi-scale Euler conjunction beat the single-scale crest set?

The shipped augmentation uses Euler solutions from one window scale (10 px = 1 km box).  A single
scale is the classic weakness of Euler deconvolution: too small a window resolves shallow contacts
and is noise-dominated, too large a window averages over several structures.  The standard remedy
(Reid et al. 1990's own discussion; Mushayandebvu et al. 2001 multiple-window schemes) is to require
a solution to be corroborated at *several* window scales.

This script builds the cross-family depth-cluster field at four window scales (10, 15, 20, 30 px),
forms five conjunction rules, and measures each resulting emission set against the incumbent with
the same blocked instrument used to ship the single-scale artifact:

  * ``N_new``  -- pixels not already emitted by the incumbent, > 200 m from the catalogue
  * ``dTP_w``  -- extra weighted credit those pixels earn on the LM instrument
  * ``eff``    -- credit per added pixel, against the metric's break-even tau(0.2778) = 0.0588
  * ``LM-cal`` -- the blocked, prevalence-calibrated instrument score of incumbent | candidate

A rule only counts as an improvement if it raises LM-cal *and* stays positive in >= 3 of 4 blocked
folds *and* keeps the emission inside the instrument's stated validity domain (< 120 k px).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import rasterio
from scipy.ndimage import distance_transform_edt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gems40.instrument import lm_score, load_live_mirror
from gems40.stack import GRAV_LAYER, MAG_LAYERS, crest, crossfamily, load_fields
from gems40.metric import ALPHA

TAU_0278 = 0.0588


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--work", default="work")
    ap.add_argument("--prior", default="/tmp/data/prior")
    ap.add_argument("--incumbent", default="top_02778_h33b2_zeros.tif")
    ap.add_argument("--scales", type=int, nargs="+", default=[10, 15, 20, 30])
    ap.add_argument("--out", default="data/evidence/multiscale.json")
    args = ap.parse_args()

    ctx = load_live_mirror(args.data)
    foot = ctx.foot
    labels = ctx.labels
    d_cat = distance_transform_edt(~labels)

    with rasterio.open(Path(args.prior) / args.incumbent) as ds:
        inc = np.isfinite(ds.read(1)) & (ds.read(1) > 0) & foot

    inc_score = lm_score(inc, ctx)
    print(f"incumbent: {int(inc.sum())} px, LM-cal {inc_score['lm_calibrated_mean']:.5f}")

    fields = {w: crossfamily(load_fields(args.work, window=w,
                                         needed=tuple(MAG_LAYERS) + (GRAV_LAYER,)), foot)
              for w in args.scales}
    crests = {w: crest(f, foot) for w, f in fields.items()}

    def robust_norm(f):                      # rank-normalise so scales are comparable
        v = np.nan_to_num(f, nan=0.0)[foot]
        order = np.argsort(np.argsort(v)) / (v.size - 1)
        out = np.full(f.shape, np.nan)
        out[foot] = order
        return out

    ranks = {w: robust_norm(f) for w, f in fields.items()}
    lo, hi = args.scales[0], args.scales[-1]

    candidates = {
        f"w{lo}-crest (shipped rule)": crests[lo],
        f"w{hi}-crest": crests[hi],
        f"w{lo} AND w{hi} crests": crests[lo] & crests[hi],
        f"all-scale crest intersection": np.logical_and.reduce([crests[w] for w in args.scales]),
        f"rank-min over scales, crest": crest(np.nanmin(np.stack([ranks[w] for w in args.scales]),
                                                        axis=0), foot),
        f"rank-mean over scales, crest": crest(np.nanmean(np.stack([ranks[w] for w in args.scales]),
                                                          axis=0), foot),
        f"w{lo} crest AND w{hi} rank>=0.90": crests[lo] &
            (np.nan_to_num(ranks[hi], nan=0.0) >= 0.90),
    }

    rows = []
    D0 = 0.2778          # the incumbent's owner-reported live score
    g_live = 12_226.0    # live-calibrated |G| (see README "Honest verdict")

    def report(name, mask):
        add = mask & ~inc & (d_cat > 2)
        n_add = int(add.sum())
        if n_add == 0:
            return
        merged = inc | add
        sc = lm_score(merged, ctx)
        d_tp = sum(d["tp_w"] for d in sc["fold_detail"].values())
        d_tp0 = sum(d["tp_w"] for d in inc_score["fold_detail"].values())
        eff = (d_tp - d_tp0) / n_add
        folds = sc["lm_calibrated_per_fold"]
        better = sum(1 for k in folds if folds[k] > inc_score["lm_calibrated_per_fold"][k])
        # projected live gain: dDTI = dn*(e - alpha*D0)*D0/T, T = D0*(beta*|G| + alpha*P)
        p_emit = float(merged.sum())
        t_live = D0 * (0.8 * g_live + 0.2 * p_emit)
        proj = n_add * (eff - ALPHA * D0) * D0 / t_live
        rows.append(dict(
            rule=name, added_px=n_add, emitted_px=int(p_emit),
            d_tp_w=round(d_tp - d_tp0, 2), efficiency=round(eff, 6),
            ratio_to_break_even=round(eff / TAU_0278, 3),
            above_break_even=bool(eff > TAU_0278),
            lm_cal=round(sc["lm_calibrated_mean"], 5),
            delta_vs_incumbent=round(sc["lm_calibrated_mean"] - inc_score["lm_calibrated_mean"], 5),
            folds_better_than_incumbent=better,
            inside_validity_domain=bool(p_emit < 120_000),
            projected_live_gain=round(proj, 5),
            lm_cal_per_fold={k: round(v, 5) for k, v in folds.items()},
        ))
        print(f"  {name:40s} N+{n_add:6d} em={int(p_emit):6d} eff={eff:9.6f} ({eff / TAU_0278:4.2f}x) "
              f"LM-cal={sc['lm_calibrated_mean']:.5f} ({rows[-1]['delta_vs_incumbent']:+.5f}, "
              f"{better}/4) proj={proj:+.5f}")

    for name, mask in candidates.items():
        report(name, mask)

    rmin = np.nanmin(np.stack([ranks[w] for w in args.scales]), axis=0)
    rmin = np.where(foot, np.nan_to_num(rmin, nan=0.0), np.nan)
    cand = crest(rmin, foot) & ~inc & (d_cat > 2)
    vals = rmin[cand]
    order = np.argsort(vals)[::-1]
    rr, cc = np.nonzero(cand)
    # ---- is the multi-scale set better than, or better *with*, the single-scale crest set? -----
    w10_add = crests[lo] & ~inc & (d_cat > 2)
    rmin_60 = np.zeros_like(inc)
    rmin_60[rr[order[:60_000]], cc[order[:60_000]]] = True
    report("rank-min crest (top 60k) UNION w10 crest", rmin_60 | w10_add)
    report("rank-min crest (top 60k) INTERSECT w10 crest", rmin_60 & w10_add)

    # ---- budget sweep of the multi-scale rank-min rule: how much of it can be shipped? ---------
    for budget in (10_000, 20_000, 30_000, 40_000, 60_000, 80_000):
        if budget >= order.size:
            continue
        keep = order[:budget]
        m = np.zeros_like(inc)
        m[rr[keep], cc[keep]] = True
        report(f"rank-min crest, top {budget // 1000}k px", m)

    rows.sort(key=lambda r: -(r["projected_live_gain"] if r["inside_validity_domain"]
                             and r["folds_better_than_incumbent"] == 4 else -1e9))
    best = rows[0]
    print(f"\nbest rule: {best['rule']}  LM-cal {best['lm_cal']} "
          f"({best['delta_vs_incumbent']:+.5f}), {best['folds_better_than_incumbent']}/4 folds, "
          f"emitted {best['emitted_px']} px")
    out = Path(args.out)
    out.write_text(json.dumps(dict(
        scales=args.scales, incumbent=args.incumbent,
        incumbent_lm_cal=round(inc_score["lm_calibrated_mean"], 5),
        incumbent_px=int(inc.sum()), tau_0278=TAU_0278,
        validity_domain_note="LM over-rewards density above ~120k emitted px (see "
                             "data/evidence/validation.json); rows above that are flagged.",
        rows=rows), indent=2))
    print(f"[out] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
