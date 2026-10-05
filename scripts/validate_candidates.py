#!/usr/bin/env python
"""Score candidate emissions with the blocked instruments, and re-verify the instrument itself.

Stage 1 (mandatory): reproduce the *known live orderings* with the LM instrument.  If the
instrument does not reproduce them, no candidate decision may be taken from it.

Stage 2: score the Euler/KDE candidates (raw, and with the catalogue-flank exclusion rule) and
print a ranked table.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import rasterio

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gems40.instrument import (cat_hidden_score, lm_score, load_cat_hidden, load_live_mirror,
                               read_binary)

# (label, path, owner-reported live DTI) -- the group's own published score ledger
LIVE_CHAIN = [
    ("H33-2-B2 (site best, 0.2778)", "top_02778_h33b2_zeros.tif", 0.2778),
    ("H27-4-R1-SOLO (0.2708)", "base_02708_h274r1_d28_nan.tif", 0.2708),
    ("D2.8 (0.2600)", "prior_d28_h195_nan.tif", 0.2600),
    ("D1.5 (0.2477)", "prior_d15_h195_nan.tif", 0.2477),
    ("T-V2-ON-D1.5 (0.2449)", "prior_tv2_d15_nan.tif", 0.2449),
    ("H19-5 solid (0.1922)", "prior_h195_solid_nan.tif", 0.1922),
    ("H16-1 (0.1911)", "prior_h16_1_nan.tif", 0.1911),
    ("H28-dotted-ridge (0.1839)", "prior_h28dotted_nan.tif", 0.1839),
    ("LATTICE-S5 (0.0904)", "prior_lattice_nan.tif", 0.0904),
    ("PLACEHOLDER (0.0445)", "prior_placeh_nan.tif", 0.0445),
]

ORDERINGS = [
    ("H33-2-B2 (site best, 0.2778)", "H27-4-R1-SOLO (0.2708)", 0.2778, 0.2708),
    ("H27-4-R1-SOLO (0.2708)", "D2.8 (0.2600)", 0.2708, 0.2600),
    ("D2.8 (0.2600)", "D1.5 (0.2477)", 0.2600, 0.2477),
    ("D1.5 (0.2477)", "T-V2-ON-D1.5 (0.2449)", 0.2477, 0.2449),
    ("H19-5 solid (0.1922)", "H16-1 (0.1911)", 0.1922, 0.1911),
    ("LATTICE-S5 (0.0904)", "PLACEHOLDER (0.0445)", 0.0904, 0.0445),
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--prior", default="/tmp/data/prior")
    ap.add_argument("--work", default="work")
    ap.add_argument("--out", default="data/evidence/validation.json")
    args = ap.parse_args()

    ctx = load_live_mirror(args.data)
    cat = load_cat_hidden(args.data)
    print(f"[LM] 4 blocked cells, flat truth px = {sum(int(c.truth.sum()) for c in ctx.cells)} "
          f"(off-catalogue SGMC)")

    # ---- Stage 1: instrument verification against the published live record --------------------
    prior_dir = Path(args.prior)
    scores: dict[str, dict] = {}
    for label, fname, lb in LIVE_CHAIN:
        p = prior_dir / fname
        if not p.exists():
            print(f"  ! missing {p}")
            continue
        m = read_binary(p)
        r = lm_score(m, ctx)
        r["leaderboard_dti"] = lb
        r["cat_hidden"] = cat_hidden_score(m, cat)["cat_hidden_mean"]
        scores[label] = r
    print(f"\n{'reference (owner-reported live)':34s} {'live':>7s} {'LM':>8s} {'LM-cal':>8s} {'N_px':>7s}")
    for label, _, lb in LIVE_CHAIN:
        if label in scores:
            r = scores[label]
            print(f"{label:34s} {lb:7.4f} {r['lm_mean']:8.5f} {r['lm_calibrated_mean']:8.5f} "
                  f"{r['emitted_pixels']:7d}")
    checks = []
    for hi, lo, lb_hi, lb_lo in ORDERINGS:
        if hi not in scores or lo not in scores:
            continue
        checks.append(dict(ordering=f"{hi} > {lo}", live=f"{lb_hi} > {lb_lo}",
                           lm=f"{scores[hi]['lm_mean']:.6f} vs {scores[lo]['lm_mean']:.6f}",
                           reproduced=bool(scores[hi]["lm_mean"] > scores[lo]["lm_mean"]),
                           lm_margin=scores[hi]["lm_mean"] - scores[lo]["lm_mean"],
                           lm_calibrated=f"{scores[hi]['lm_calibrated_mean']:.6f} vs "
                                         f"{scores[lo]['lm_calibrated_mean']:.6f}",
                           reproduced_calibrated=bool(scores[hi]["lm_calibrated_mean"] >
                                                      scores[lo]["lm_calibrated_mean"])))
    n_ok = sum(c["reproduced_calibrated"] for c in checks)
    print(f"\n[instrument] live orderings reproduced by LM-calibrated: {n_ok}/{len(checks)}")
    for c in checks:
        print(f"   {'OK ' if c['reproduced_calibrated'] else 'MISS'} {c['ordering']:46s} "
              f"LM {c['lm_calibrated']}")

    report = dict(instrument=dict(
        name="LM (spatially blocked, off-catalogue SGMC truth, prevalence-calibrated)",
        n_orderings=len(checks), n_reproduced=n_ok, orderings=checks,
        references={k: dict(live=v["leaderboard_dti"], lm=v["lm_mean"],
                            lm_calibrated=v["lm_calibrated_mean"],
                            cat_hidden=v["cat_hidden"], emitted_pixels=v["emitted_pixels"])
                    for k, v in scores.items()}),
        candidates={})

    # ---- Stage 2: Euler/KDE candidates ---------------------------------------------------------
    work = Path(args.work)
    fields = {}
    for npz in sorted(work.glob("euler_*.npz")):
        d = np.load(npz, allow_pickle=True)
        fields[str(d["layer"])] = (npz.name, {k: d[k] for k in d.files if k != "field"}, d["field"])

    print("\n[candidates] %d Euler/KDE fields: %s" % (len(fields), ", ".join(fields)))

    def support_from_field(field: np.ndarray, foot_ctx, n_target: int) -> np.ndarray:
        flat = field[foot_ctx]
        if n_target >= flat.size:
            return np.isfinite(field) & (field > 0)
        thr = np.partition(flat, flat.size - n_target)[flat.size - n_target]
        return np.nan_to_num(field, nan=-1.0) >= thr

    foot = ctx.foot
    results = []
    for layer, (name, meta, field) in fields.items():
        for n_target in (30_000, 44_090, 60_069, 90_000):
            m = support_from_field(field, foot, n_target)
            for prune in (None, 2):
                mp = m
                if prune is not None:
                    from scipy.ndimage import distance_transform_edt
                    d_cat = distance_transform_edt(~ctx.labels)
                    mp = m & (d_cat > prune)
                r = lm_score(mp, ctx)
                rc = cat_hidden_score(mp, cat)
                results.append(dict(layer=layer, n_target=n_target, prune=prune,
                                    n_px=int(mp.sum()), lm=r["lm_calibrated_mean"],
                                    lm_raw=r["lm_mean"], cat_hidden=rc["cat_hidden_mean"],
                                    on_cat=r["on_catalogue_pixels"]))
    results.sort(key=lambda d: -d["lm"])
    print(f"\n{'layer':14s} {'target':>7s} {'prune':>5s} {'N_px':>7s} {'LM-cal':>8s} {'CAT-HID':>8s} {'on_cat':>6s}")
    for d in results:
        print(f"{d['layer']:14s} {d['n_target']:7d} {str(d['prune']):>5s} {d['n_px']:7d} "
              f"{d['lm']:8.5f} {d['cat_hidden']:8.5f} {d['on_cat']:6d}")
    report["candidates"] = results
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2))
    print(f"\n[out] {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
