#!/usr/bin/env python
"""Refinement test: can a crest-thinned (NMS) version of the Euler depth-cluster field beat the
break-even marginal efficiency on top of the incumbent?

Motivation (measured, not assumed): the raw KDE support is *almost disjoint* from the incumbent's
ridge family, so it is new information rather than redundancy -- but its credit-per-pixel is below
the break-even tau.  A KDE ridge is a 2-D band; the metric pays per emitted pixel, so concentrating
the same cloud onto 1-px crests (non-maximum suppression, as used for the group's other ridge
surfaces) is the natural way to raise credit-per-pixel.

All numbers are measured on the LM instrument with the same break-even test as `union_test.py`.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import distance_transform_edt, maximum_filter

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gems40.instrument import load_live_mirror, read_binary
from gems40.metric import ALPHA

sys.path.insert(0, str(Path(__file__).resolve().parent))
from union_test import tau, tp_w_matrix  # noqa: E402


def nms_crest(field: np.ndarray, size: int = 3, threshold_quantile: float = 0.0) -> np.ndarray:
    """Local-maximum crests of a continuous field (ties broken by the maximum filter)."""
    f = np.nan_to_num(field, nan=-np.inf)
    mx = maximum_filter(f, size=size, mode="nearest")
    crest = (f >= mx) & np.isfinite(field)
    if threshold_quantile > 0:
        thr = np.nanquantile(field[np.isfinite(field)], threshold_quantile)
        crest &= field >= thr
    return crest


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--prior", default="/tmp/data/prior")
    ap.add_argument("--work", default="work")
    ap.add_argument("--incumbent", default="top_02778_h33b2_zeros.tif")
    ap.add_argument("--out", default="data/evidence/nms_refine.json")
    args = ap.parse_args()

    ctx = load_live_mirror(args.data)
    d_cat = distance_transform_edt(~ctx.labels)
    inc = read_binary(Path(args.prior) / args.incumbent) & ctx.foot
    tp_inc, g_tot, n_inc = tp_w_matrix(inc, ctx)
    brk = tau(0.2778)
    print(f"incumbent N={n_inc} TP_w={tp_inc:.1f} | break-even tau(0.2778)={brk:.4f} credit/px")

    fields = {}
    for npz in sorted(Path(args.work).glob("euler_*.npz")):
        d = np.load(npz, allow_pickle=True)
        fields[str(d["layer"])] = d["field"].astype(np.float64)
    foot = ctx.foot

    # normalised family combination: every layer is scaled by its own in-footprint maximum
    def norm(f):
        m = np.nanmax(f[foot])
        return np.where(foot, f / m, np.nan)

    mags = [norm(fields[k]) for k in ("rtp", "tmi", "mag_anom") if k in fields]
    mag = np.nanmean(np.stack(mags), axis=0) if mags else None
    grav = norm(fields["iso_grav_anom"]) if "iso_grav_anom" in fields else None
    conj = np.fmin(mag, grav) if (mag is not None and grav is not None) else None
    combos = dict(mag_mean=mag, gravity=grav, conjunction=conj)

    rows = []
    for name, f in combos.items():
        if f is None:
            continue
        # (a) plain top-N threshold on the KDE
        for n_t in (5_000, 10_000, 20_000):
            flat = f[foot]
            thr = np.partition(flat, flat.size - n_t)[flat.size - n_t]
            add = (f >= thr) & foot & ~inc & (d_cat > 2)
            n_add = int(add.sum())
            if n_add == 0:
                continue
            tp_u, _, _ = tp_w_matrix(inc | add, ctx)
            eff = (tp_u - tp_inc) / n_add
            rows.append(dict(candidate=f"{name}:top{n_t}", kind="kde_top",
                             added_px=n_add, efficiency=eff, ratio=eff / brk))
        # (b) 1-px NMS crests of the KDE, then a top-N cap on crest density
        for n_t in (5_000, 10_000, 20_000, 40_000):
            flat = f[foot]
            thr = np.partition(flat, flat.size - n_t)[flat.size - n_t]
            crest = nms_crest(f, size=3) & (f >= thr) & foot
            add = crest & ~inc & (d_cat > 2)
            n_add = int(add.sum())
            if n_add == 0:
                continue
            tp_u, _, _ = tp_w_matrix(inc | add, ctx)
            eff = (tp_u - tp_inc) / n_add
            rows.append(dict(candidate=f"{name}:nms{n_t}", kind="kde_nms",
                             added_px=n_add, efficiency=eff, ratio=eff / brk))
    rows.sort(key=lambda r: -r["efficiency"])
    print(f"\n{'candidate':34s} {'kind':9s} {'added_px':>9s} {'eff':>9s} {'ratio':>6s} > tau?")
    for r in rows:
        print(f"{r['candidate']:34s} {r['kind']:9s} {r['added_px']:9d} {r['efficiency']:9.6f} "
              f"{r['ratio']:6.2f} {'YES' if r['efficiency'] > brk else 'no'}")
    Path(args.out).write_text(json.dumps(dict(incumbent=args.incumbent, incumbent_tp_w=tp_inc,
                                              incumbent_n_px=n_inc, tau_0278=brk, rows=rows),
                                         indent=2))
    print(f"\n[out] {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
