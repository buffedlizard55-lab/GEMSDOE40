#!/usr/bin/env python
"""Marginal-value test: does Euler/KDE depth-clustered mass add credit on top of the incumbent?

The live metric is

    DTI = TP_w / (alpha*(TP_w + FP_w) + beta*|G| + eps)

so, at an operating point D0, an addition of ``dn`` prediction pixels that earns ``dTP_w`` extra
weighted credit changes the score by

    d DTI > 0   <=>   dTP_w / dn  >  tau(D0) = alpha*D0 / (1 - alpha*D0)

(Sibling repository ``GEMSDOE28``, ``knowledge/43_hypotheses_heatflow_euler.md`` §1.1; the
identity is re-derived and unit-checked in ``tests/test_metric.py``.)

This script measures ``dTP_w`` and ``dn`` **on the LM instrument** (off-catalogue SGMC truth,
4 spatially blocked quadrants, prevalence calibrated) for each candidate addition, and reports
the marginal efficiency in credit-per-pixel against the break-even tau at the incumbent's own
measured score.  Nothing is shipped from this script; it produces the evidence that decides.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import distance_transform_edt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gems40.instrument import load_live_mirror, read_binary
from gems40.metric import ALPHA, BETA, EPS, kernel

BREAK_EVEN_LIVE = 0.3195   # the score to beat


def tau(d0: float) -> float:
    return ALPHA * d0 / (1.0 - ALPHA * d0)


def tp_w_matrix(mask: np.ndarray, ctx) -> tuple[float, float, int]:
    """Return (TP_w, |G|, N) on the LM cells for one mask (as used by the LM instrument)."""
    m = np.asarray(mask, bool) & ctx.foot
    tot_tp, tot_g, tot_n = 0.0, 0.0, 0
    for cell in ctx.cells:
        sl = cell.bbox
        p = m[sl] & cell.domain
        g = cell.truth
        n_p, n_g = int(p.sum()), int(g.sum())
        tot_n += n_p
        tot_g += n_g
        if n_p == 0 or n_g == 0:
            continue
        dp = distance_transform_edt(~p)
        dg = distance_transform_edt(~g)
        tp_p = float(kernel(dp[g]).sum())
        tp_g = float(kernel(dg[p]).sum())
        tot_tp += 0.5 * (tp_p + tp_g)
    return tot_tp, tot_g, tot_n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--prior", default="/tmp/data/prior")
    ap.add_argument("--work", default="work")
    ap.add_argument("--incumbent", default="top_02778_h33b2_zeros.tif")
    ap.add_argument("--out", default="data/evidence/union_test.json")
    args = ap.parse_args()

    ctx = load_live_mirror(args.data)
    d_cat = distance_transform_edt(~ctx.labels)

    inc = read_binary(Path(args.prior) / args.incumbent) & ctx.foot
    tp_inc, g_tot, n_inc = tp_w_matrix(inc, ctx)
    print(f"incumbent {args.incumbent}: N={n_inc} px, TP_w={tp_inc:.1f} (LM truth |G|={int(g_tot)})")
    print(f"break-even marginal efficiency at live 0.2778 = {tau(0.2778):.4f} credit/px; "
          f"at live 0.3195 = {tau(BREAK_EVEN_LIVE):.4f} credit/px")

    rows = []
    fields = {}
    for npz in sorted(Path(args.work).glob("euler_*.npz")):
        d = np.load(npz, allow_pickle=True)
        fields[str(d["layer"])] = d["field"]
    mag = [fields[k] for k in ("rtp", "tmi", "mag_anom") if k in fields]
    grav = fields.get("iso_grav_anom")
    foot = ctx.foot

    def support(field, n_target):
        flat = field[foot]
        thr = np.partition(flat, flat.size - n_target)[flat.size - n_target]
        return np.nan_to_num(field, nan=-1.0) >= thr

    def report(name, add_mask):
        add = add_mask & ~inc & (d_cat > 2)          # only genuinely new, off-catalogue-flank mass
        n_add = int(add.sum())
        if n_add == 0:
            return
        tp_u, _, _ = tp_w_matrix(inc | add, ctx)
        d_tp = tp_u - tp_inc
        eff = d_tp / n_add
        rows.append(dict(candidate=name, added_px=n_add, d_tp_w=round(d_tp, 2),
                         efficiency=round(eff, 6), break_even_0278=round(tau(0.2778), 6),
                         above_break_even=bool(eff > tau(0.2778)),
                         ratio_to_break_even=round(eff / tau(0.2778), 3)))

    for layer, f in fields.items():
        for n_t in (10_000, 20_000, 30_000, 44_090):
            report(f"{layer}:top{n_t}", support(f, n_t))
    if len(mag) > 0 and grav is not None:
        m = np.minimum(np.nan_to_num(mag[0], nan=-1), np.nan_to_num(grav, nan=-1))
        report("mag_x_grav_conjunction_top20000", support(m, 20_000))
        report("mag_x_grav_conjunction_top40000", support(m, 40_000))

    rows.sort(key=lambda r: -r["efficiency"])
    print(f"\n{'candidate':38s} {'added_px':>9s} {'dTP_w':>8s} {'eff':>9s} {'ratio':>6s}  > tau?")
    for r in rows:
        print(f"{r['candidate']:38s} {r['added_px']:9d} {r['d_tp_w']:8.1f} {r['efficiency']:9.6f} "
              f"{r['ratio_to_break_even']:6.2f}  {'YES' if r['above_break_even'] else 'no'}")
    Path(args.out).write_text(json.dumps(dict(
        incumbent=args.incumbent, incumbent_n_px=n_inc, incumbent_tp_w=tp_inc,
        lm_truth_px=int(g_tot), tau_live_0278=tau(0.2778), tau_live_03195=tau(BREAK_EVEN_LIVE),
        rows=rows), indent=2))
    print(f"\n[out] {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
