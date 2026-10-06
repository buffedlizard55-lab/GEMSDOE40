#!/usr/bin/env python
"""Final audit of the shipped artifacts: validation on the blocked instruments + novelty check.

Novelty rule (fixed before looking at the numbers, so it cannot be tuned after the fact):

    near-duplicate  <=>  Jaccard(support, prior support) > 0.50
                     OR  |Pearson r| of the two rasters (footprint-restricted, mass-normalised)
                         > 0.80

Both thresholds are stated in the brief's terms ("hash and correlate it against every prior
submission's raw output and refuse to call it new if it's a near-duplicate").  The audit prints
the full matrix, not just the verdict, so any reader can re-judge with different thresholds.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import rasterio

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gems40.instrument import cat_hidden_score, lm_score, load_cat_hidden, load_live_mirror, read_binary

JACCARD_DUP = 0.50
PEARSON_DUP = 0.80


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_full(path: Path, foot: np.ndarray) -> np.ndarray:
    with rasterio.open(path) as ds:
        a = ds.read(1).astype(np.float64)
    return np.where(np.isfinite(a) & foot, a, np.nan)


def pearson(a: np.ndarray, b: np.ndarray, foot: np.ndarray) -> float:
    x = np.nan_to_num(a, nan=0.0)[foot]
    y = np.nan_to_num(b, nan=0.0)[foot]
    if x.std() == 0 or y.std() == 0:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def spear(a: np.ndarray, b: np.ndarray, foot: np.ndarray, stride: int = 7) -> float:
    x = np.nan_to_num(a, nan=0.0)[foot][::stride]
    y = np.nan_to_num(b, nan=0.0)[foot][::stride]
    rx = np.argsort(np.argsort(x)).astype(np.float64)
    ry = np.argsort(np.argsort(y)).astype(np.float64)
    if rx.std() == 0 or ry.std() == 0:
        return 0.0
    return float(np.corrcoef(rx, ry)[0, 1])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--downloads", default="docs/downloads")
    ap.add_argument("--prior", default="/tmp/data/prior")
    ap.add_argument("--data", default="data")
    ap.add_argument("--outdir", default="data/evidence")
    args = ap.parse_args()

    ctx = load_live_mirror(args.data)
    cat = load_cat_hidden(args.data)
    foot = ctx.foot
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ddir = Path(args.downloads)
    artifacts = sorted(p for p in ddir.glob("*.tif"))
    priors = sorted(p for p in Path(args.prior).glob("*.tif") if p.name != "sgmc_faults_u8.tif")

    # ---- hashes ------------------------------------------------------------------------------
    hashes = dict(artifacts={p.name: sha256(p) for p in artifacts},
                  priors={p.name: sha256(p) for p in priors})
    dup_hash = {p.name: [q.name for q in priors if hashes["priors"][q.name] == hashes["artifacts"][p.name]]
                for p in artifacts}

    # ---- validation --------------------------------------------------------------------------
    inc = read_binary(Path(args.prior) / "top_02778_h33b2_zeros.tif") & foot
    base = lm_score(inc, ctx)
    base_ch = cat_hidden_score(inc, cat)
    val = dict(incumbent=dict(name="top_02778_h33b2_zeros.tif (owner-reported live 0.2778)",
                              lm_calibrated=base["lm_calibrated_mean"],
                              lm_raw=base["lm_mean"], cat_hidden=base_ch["cat_hidden_mean"],
                              n_px=base["emitted_pixels"],
                              folds=base["lm_calibrated_per_fold"]),
               artifacts={})
    for p in artifacts:
        m = read_binary(p) & foot
        r = lm_score(m, ctx)
        ch = cat_hidden_score(m, cat)
        folds = r["lm_calibrated_per_fold"]
        bf = base["lm_calibrated_per_fold"]
        val["artifacts"][p.name] = dict(
            lm_calibrated=r["lm_calibrated_mean"], lm_raw=r["lm_mean"],
            cat_hidden=ch["cat_hidden_mean"], n_px=r["emitted_pixels"],
            on_catalogue_px=r["on_catalogue_pixels"], folds=folds,
            delta_vs_incumbent=r["lm_calibrated_mean"] - base["lm_calibrated_mean"],
            folds_better_than_incumbent=int(sum(folds[k] > bf[k] for k in folds)),
            beats_incumbent=bool(r["lm_calibrated_mean"] > base["lm_calibrated_mean"]),
            all_folds_better=bool(all(folds[k] > bf[k] for k in folds)))

    # ---- novelty -----------------------------------------------------------------------------
    prior_fields = {p.name: read_full(p, foot) for p in priors}
    rows = []
    for p in artifacts:
        f = read_full(p, foot)
        for q_name, qf in prior_fields.items():
            inter = float(((f > 0) & (qf > 0)).sum())
            union = float(((f > 0) | (qf > 0)).sum())
            jac = inter / union if union else 0.0
            pr = pearson(f, qf, foot)
            sp = spear(f, qf, foot)
            rows.append(dict(artifact=p.name, prior=q_name, jaccard=round(jac, 6),
                             pearson=round(pr, 6), spearman=round(sp, 6),
                             near_duplicate=bool(jac > JACCARD_DUP or abs(pr) > PEARSON_DUP)))
    summary = {}
    for p in artifacts:
        sub = [r for r in rows if r["artifact"] == p.name]
        summary[p.name] = dict(
            max_jaccard=max(r["jaccard"] for r in sub),
            max_abs_pearson=max(abs(r["pearson"]) for r in sub),
            n_priors_compared=len(sub),
            n_near_duplicates=sum(r["near_duplicate"] for r in sub),
            verdict=("NEAR-DUPLICATE" if any(r["near_duplicate"] for r in sub) else "NEW"))
    novelty = dict(rule=dict(jaccard_duplicate_threshold=JACCARD_DUP,
                             abs_pearson_duplicate_threshold=PEARSON_DUP,
                             note="verdict NEW requires both: Jaccard <= 0.5 and |Pearson| <= 0.8 "
                                  "against every prior artifact"),
                   hashes=hashes, identical_to_prior=dup_hash, matrix=rows, summary=summary)

    (outdir / "shipped_validation.json").write_text(json.dumps(val, indent=2))
    (outdir / "novelty.json").write_text(json.dumps(novelty, indent=2))

    print(f"{'artifact':64s} {'N':>7s} {'LM-cal':>8s} {'incumbent':>10s} {'folds>inc':>9s}")
    for name, d in val["artifacts"].items():
        print(f"{name:64s} {d['n_px']:7d} {d['lm_calibrated']:8.5f} "
              f"{d['delta_vs_incumbent']:+10.5f} {d['folds_better_than_incumbent']:9d}")
    print(f"\n{'incumbent reference':64s} {base['emitted_pixels']:7d} "
          f"{base['lm_calibrated_mean']:8.5f}")
    print(f"\n{'artifact':64s} {'maxJacc':>8s} {'max|r|':>7s} {'dups':>5s}  verdict")
    for name, d in summary.items():
        print(f"{name:64s} {d['max_jaccard']:8.4f} {d['max_abs_pearson']:7.4f} "
              f"{d['n_near_duplicates']:5d}  {d['verdict']}")
    print(f"\n[out] {outdir/'shipped_validation.json'} , {outdir/'novelty.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
