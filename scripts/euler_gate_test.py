#!/usr/bin/env python
"""Test: can the Euler depth-cluster field *gate* (prune) the incumbent emission?

Physically motivated rule under test -- the "concealed-contact licence":

    keep an incumbent dot only if a shallow (SI = 0) Euler contact cluster is present within
    r pixels of it

Rationale: lineaments without any shallow potential-field contact underneath are the ones most
likely to be non-fault artefacts (roads, strandlines, unit boundaries), which is the class of
false positive the group has repeatedly failed to remove.  Removal is charged by the metric at
alpha per pixel, so the gain condition for a removal of dn dots that costs dS weighted credit is
the same break-even test as for an addition (Sibling repo GEMSDOE28 knowledge/43 §2.2).

Everything is measured on the LM instrument.  A prune is only "promoted" if it improves
LM-calibrated mean in **all four** blocked folds.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import binary_dilation, distance_transform_edt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gems40.instrument import lm_score, load_live_mirror, read_binary


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--prior", default="/tmp/data/prior")
    ap.add_argument("--work", default="work")
    ap.add_argument("--incumbent", default="top_02778_h33b2_zeros.tif")
    ap.add_argument("--out", default="data/evidence/euler_gate.json")
    args = ap.parse_args()

    ctx = load_live_mirror(args.data)
    inc = read_binary(Path(args.prior) / args.incumbent) & ctx.foot
    base = lm_score(inc, ctx)
    print(f"incumbent: N={base['emitted_pixels']} LM-cal={base['lm_calibrated_mean']:.6f} "
          f"folds={ {k: round(v, 4) for k, v in base['lm_calibrated_per_fold'].items()} }")

    fields = {}
    for npz in sorted(Path(args.work).glob("euler_*.npz")):
        d = np.load(npz, allow_pickle=True)
        fields[str(d["layer"])] = d["field"]
    foot = ctx.foot

    mag_keys = [k for k in ("rtp", "tmi", "mag_anom") if k in fields]
    if mag_keys:
        m = np.minimum.reduce([np.nan_to_num(fields[k], nan=-1.0) for k in mag_keys])
        fields["mag_union"] = np.where(foot, m, np.nan)

    rows = []
    for layer, f in fields.items():
        for n_sup in (20_000, 44_090, 90_000):
            flat = f[foot]
            thr = np.partition(flat, flat.size - n_sup)[flat.size - n_sup]
            sup = (np.nan_to_num(f, nan=-1.0) >= thr) & foot
            for r in (1, 2, 3):
                gate = binary_dilation(sup, iterations=r) & foot
                keep = inc & gate
                n_keep = int(keep.sum())
                if n_keep < 5_000:
                    continue
                r_ = lm_score(keep, ctx)
                folds = r_["lm_calibrated_per_fold"]
                base_folds = base["lm_calibrated_per_fold"]
                all_better = all(folds[k] > base_folds[k] for k in folds)
                rows.append(dict(layer=layer, euler_support_px=n_sup, dilate_px=r,
                                 kept_px=n_keep, removed_px=int(inc.sum()) - n_keep,
                                 lm_cal=r_["lm_calibrated_mean"],
                                 delta=r_["lm_calibrated_mean"] - base["lm_calibrated_mean"],
                                 folds_positive=int(sum(folds[k] > base_folds[k] for k in folds)),
                                 all_folds_better=bool(all_better)))
    rows.sort(key=lambda d: -d["delta"])
    print(f"\n{'euler field':12s} {'support':>8s} {'dil':>4s} {'kept':>7s} {'removed':>8s} "
          f"{'LM-cal':>8s} {'delta':>9s} {'folds+':>6s}")
    for d in rows[:25]:
        print(f"{d['layer']:12s} {d['euler_support_px']:8d} {d['dilate_px']:4d} {d['kept_px']:7d} "
              f"{d['removed_px']:8d} {d['lm_cal']:8.5f} {d['delta']:+9.5f} {d['folds_positive']:6d}")
    Path(args.out).write_text(json.dumps(dict(
        incumbent=args.incumbent, incumbent_lm_cal=base["lm_calibrated_mean"],
        incumbent_n_px=base["emitted_pixels"], rows=rows), indent=2))
    print(f"\n[out] {args.out}")
    best = rows[0] if rows else None
    if best and best["all_folds_better"]:
        print(f"PROMOTABLE: {best}")
    else:
        print("No prune improved LM-calibrated mean in all four blocked folds.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
