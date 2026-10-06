#!/usr/bin/env python3
"""Empirical study of the prior-submission corpus on the competition grid.

Answers, with measurements rather than assertions:

1. What is the *geometry* of each scored prior (emitted pixels, dot lattice
   spacing, redundancy inside the 300 m kernel)?
2. How does each prior score on the two surrogate instruments, both at its own
   emitted mass and at matched budgets (so the comparison isolates placement
   from size)?
3. What do the controls score (uniform random dots of the same mass; the
   catalogue flank ring)?
4. How strongly do the surrogate scores and the emitted mass correlate with the
   official (owner-reported) leaderboard scores?

Output: ``docs/data/prior_corpus_study.json`` plus a console table.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe40.corpus import PRIORS, load_field, local_path, sha256  # noqa: E402
from gemsdoe40.instrument import (  # noqa: E402
    binary_dti_kdtree, flank_ring, isolated_catalogue_components, random_dots,
    sgmc_off_catalogue,
)
from gemsdoe40.lattice import redundancy_stats  # noqa: E402

BUDGETS = (10_000, 20_000, 30_000, 40_000, 50_000, 60_000, 80_000, 120_000)
SGMC_PATH = ROOT / "data" / "external" / "sgmc" / "derived_sgmc_faults_100m_u8.tif"


def prefix_curve(values: np.ndarray, truth: np.ndarray, mask: np.ndarray, valid: np.ndarray,
                 budgets) -> list[dict]:
    """Exact masked DTI of the top-N pixels by value, for every budget in ``budgets``."""
    flat = np.asarray(values, dtype=np.float64).ravel()
    order = np.argsort(-flat, kind="stable")
    n_pos = int((flat > 0).sum())
    out = []
    for n in budgets:
        k = int(min(n, n_pos))
        dots = np.zeros(flat.size, dtype=bool)
        if k > 0:
            dots[order[:k]] = True
        dots = dots.reshape(values.shape)
        r = binary_dti_kdtree(dots, truth, mask, valid)
        out.append(dict(n_px=int(k), **{kk: r[kk] for kk in
                                        ("tp", "fp", "dti", "mean_credit", "kappa")}))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", default=str(ROOT / "data" / "sample_submission.tif"))
    ap.add_argument("--labels", default=str(ROOT / "data" / "labels.tif"))
    ap.add_argument("--sgmc", default=str(SGMC_PATH))
    ap.add_argument("--out", default=str(ROOT / "docs" / "data" / "prior_corpus_study.json"))
    ap.add_argument("--budgets", default=",".join(str(b) for b in BUDGETS))
    args = ap.parse_args()
    budgets = tuple(int(b) for b in args.budgets.split(","))

    t0 = time.time()
    with rasterio.open(args.sample) as s:
        valid = np.isfinite(s.read(1))
    with rasterio.open(args.labels) as s:
        catalogue = s.read(1) == 1
    with rasterio.open(args.sgmc) as s:
        sgmc = s.read(1)
    sgmc_off = sgmc_off_catalogue(sgmc, catalogue, valid)
    iso = isolated_catalogue_components(catalogue, valid, min_separation_px=6, min_size_px=12)
    iso3 = isolated_catalogue_components(catalogue, valid, min_separation_px=3, min_size_px=12)
    print(f"footprint {int(valid.sum()):,}  catalogue {int(catalogue.sum()):,}  "
          f"sgmc-off-catalogue {int(sgmc_off.sum()):,}  "
          f"isolated(>=6px) {iso['truth_px']:,} px in {iso['n_components']} components", flush=True)

    instruments = {
        "sgmc_off_catalogue": dict(truth=sgmc_off, mask=catalogue),
        "isolated_components": dict(truth=iso["truth"], mask=iso["mask"]),
        "isolated_components_sep3": dict(truth=iso3["truth"], mask=iso3["mask"]),
    }

    report = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "grid": dict(shape=list(valid.shape), footprint_px=int(valid.sum()),
                     catalogue_px=int(catalogue.sum())),
        "surrogates": {
            "sgmc_off_catalogue": dict(truth_px=int(sgmc_off.sum()),
                                       description="USGS SGMC faults absent from the catalogue"),
            "isolated_components": dict(truth_px=int(iso["truth_px"]),
                                        min_separation_px=6, min_size_px=12,
                                        description="catalogue components >=6 px from all others"),
            "isolated_components_sep3": dict(truth_px=int(iso3["truth_px"]),
                                             min_separation_px=3, min_size_px=12,
                                             description="same, 3 px separation (sensitivity)"),
        },
        "budgets_px": list(budgets),
        "priors": [],
        "controls": {},
    }

    for p in PRIORS:
        path = local_path(p)
        if not path.exists():
            print(f"MISSING {p.key} ({path}) - run scripts/fetch_prior_corpus.sh", flush=True)
            continue
        field = load_field(path, valid)
        n_px = int((field > 0).sum())
        geom = redundancy_stats(field > 0)
        rows = dict(key=p.key, repo=p.repo, path=p.path, score_owner_reported=p.score,
                    sha256_file=sha256(path), emitted_px=n_px, geometry=geom,
                    use_for_calibration=p.use_for_calibration,
                    instruments={})
        for iname, cfg in instruments.items():
            full = binary_dti_kdtree(field > 0, cfg["truth"], cfg["mask"], valid)
            curve = prefix_curve(field, cfg["truth"], cfg["mask"], valid, budgets)
            rows["instruments"][iname] = dict(full_mass=full, curve=curve)
        report["priors"].append(rows)
        s1 = rows["instruments"]["sgmc_off_catalogue"]
        s3 = rows["instruments"]["isolated_components"]
        at40 = [c["dti"] for c in s1["curve"] if c["n_px"] == 40000]
        at40_txt = f"{at40[0]:.4f}" if at40 else "n/a"
        print(f"{p.key:32s} score={p.score:.4f} px={n_px:7d} "
              f"sgmc_full={s1['full_mass']['dti']:.4f} sgmc@40k={at40_txt} "
              f"iso_full={s3['full_mass']['dti']:.4f} ({time.time() - t0:.0f}s)", flush=True)

    # controls -------------------------------------------------------------
    for iname, cfg in instruments.items():
        ctrl = {}
        for n in (20_000, 40_000):
            rs, rk = [], []
            for seed in range(3):
                dots = random_dots(valid, cfg["mask"], n, seed=seed)
                r = binary_dti_kdtree(dots, cfg["truth"], cfg["mask"], valid)
                rs.append(r["dti"]); rk.append(r["kappa"])
            ctrl[f"random_{n}"] = dict(dti_mean=float(np.mean(rs)), dti_max=float(np.max(rs)),
                                       dti_min=float(np.min(rs)), kappa_mean=float(np.mean(rk)))
        ring = flank_ring(catalogue, valid, inner_px=0, outer_px=3, mask_outside=True)
        ring &= ~catalogue
        ring_vals = np.where(ring, 1.0, 0.0)
        ctrl["catalogue_flank_ring_0_3px"] = dict(
            n_px=int(ring.sum()),
            **{k: v for k, v in binary_dti_kdtree(ring, cfg["truth"], cfg["mask"], valid).items()},
            curve=prefix_curve(ring_vals, cfg["truth"], cfg["mask"], valid, budgets),
        )
        report["controls"][iname] = ctrl
        print(f"[control {iname}] random@40k dti={ctrl['random_40000']['dti_mean']:.4f} "
              f"kappa={ctrl['random_40000']['kappa_mean']:.4f}  flank ring px={int(ring.sum())} "
              f"dti={ctrl['catalogue_flank_ring_0_3px']['dti']:.4f}", flush=True)

    # correlations ---------------------------------------------------------
    try:
        from scipy.stats import spearmanr
        scores = np.array([r["score_owner_reported"] for r in report["priors"]], dtype=float)
        corr = {}
        for label, fn in (
            ("emitted_px", lambda r: r["emitted_px"]),
            ("sgmc_full_dti", lambda r: r["instruments"]["sgmc_off_catalogue"]["full_mass"]["dti"]),
            ("sgmc_curve_mean_dti", lambda r: float(np.mean([c["dti"] for c in
                                                             r["instruments"]["sgmc_off_catalogue"]["curve"]]))),
            ("sgmc_at_equivalent_budget", None),
        ):
            if fn is None:
                continue
            x = np.array([fn(r) for r in report["priors"]], dtype=float)
            rho = spearmanr(x, scores)
            corr[label] = dict(rho=float(rho.statistic), p=float(rho.pvalue), n=len(x))
        report["correlations_vs_official"] = corr
        for k, v in corr.items():
            print(f"spearman {k:26s} rho={v['rho']:+.3f} p={v['p']:.3f} n={v['n']}")
    except Exception as exc:  # pragma: no cover
        report["correlations_vs_official"] = {"error": str(exc)}

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=1))
    print(f"wrote {out}  ({time.time() - t0:.0f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
