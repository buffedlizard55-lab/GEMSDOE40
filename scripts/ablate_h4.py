#!/usr/bin/env python3
"""Ablation of the H4 ranking configuration against the two surrogate instruments.

Every configuration is assembled by the *production* code path
(:func:`gemsdoe40.euler_h4_line.assemble_ranking`) and scored by the production lattice
and instruments, so the table below is exactly what the pipeline would ship.
The calibration fitted to the prior corpus maps each prefix curve to a modelled
distance-weighted Tversky index; the surrogate columns are measured directly.

Writes ``docs/data/h4_ablation.json``.
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

from gemsdoe40 import BANDS  # noqa: E402
from gemsdoe40.calibration import fit_model, predict_score  # noqa: E402
from gemsdoe40.euler import EulerCloud  # noqa: E402
from gemsdoe40.euler_h4_line import (  # noqa: E402
    H4Config, assemble_ranking, curvature_evidence, horizontal_gradient_magnitude, tilt_angle,
)
from gemsdoe40.grid import footprint_from_sample, read_band, read_labels  # noqa: E402
from gemsdoe40.instrument import (  # noqa: E402
    binary_dti_kdtree, isolated_catalogue_components, sgmc_off_catalogue,
)
from gemsdoe40.lattice import greedy_min_separation  # noqa: E402

WORK = ROOT / "work" / "h4"

CONFIGS: dict[str, dict] = {
    # name: H4Config overrides
    "flank26_lineKDE": dict(use_flank_band=True, flank_lo_px=2.0, flank_hi_px=6.0,
                            flank_mix=0.25, flank_weight=0.75, concord_in_line=False),
    "flank26_linePRIMARY": dict(use_flank_band=True, flank_lo_px=2.0, flank_hi_px=6.0,
                                flank_mix=0.25, flank_weight=0.75, concord_in_line=True),
    "flank36_lineKDE": dict(use_flank_band=True, flank_lo_px=3.0, flank_hi_px=6.0,
                            flank_mix=0.25, flank_weight=0.75, concord_in_line=False),
    "flank26_lineKDE_mix10": dict(use_flank_band=True, flank_lo_px=2.0, flank_hi_px=6.0,
                                  flank_mix=0.10, flank_weight=0.90, concord_in_line=False),
    "flank26_lineKDE_half1": dict(use_flank_band=True, flank_lo_px=2.0, flank_hi_px=6.0,
                                  flank_mix=0.25, flank_weight=0.75, concord_in_line=False,
                                  line_half_len=1),
    "flank26_lineKDE_half6": dict(use_flank_band=True, flank_lo_px=2.0, flank_hi_px=6.0,
                                  flank_mix=0.25, flank_weight=0.75, concord_in_line=False,
                                  line_half_len=6),
    "flank26_lineKDE_wcoh0": dict(use_flank_band=True, flank_lo_px=2.0, flank_hi_px=6.0,
                                  flank_mix=0.25, flank_weight=0.75, concord_in_line=False,
                                  w_coherence=0.0),
    "flank28_lineKDE": dict(use_flank_band=True, flank_lo_px=2.0, flank_hi_px=8.0,
                            flank_mix=0.25, flank_weight=0.75, concord_in_line=False),
    "flank26_lineKDE_sig32": dict(use_flank_band=True, flank_lo_px=2.0, flank_hi_px=6.0,
                                  flank_mix=0.25, flank_weight=0.75, concord_in_line=False,
                                  sigma_px=3.2),
    "band_only_mix0": dict(use_flank_band=True, flank_lo_px=2.0, flank_hi_px=6.0,
                           flank_mix=0.0, flank_weight=1.0, concord_in_line=False),
}


def load_cloud(path: Path) -> EulerCloud:
    z = np.load(path)
    return EulerCloud(path.stem, 0.0, -1, -1, z["p_row"], z["p_col"], z["p_depth_m"],
                      z["p_rel_se"], np.zeros_like(z["p_row"]), z["p_analytic"], {})


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--budgets", default="10000,20000,30000,35000,40000,50000,60000")
    ap.add_argument("--spacings", default="2.8284,4.0")
    ap.add_argument("--out", default=str(ROOT / "docs" / "data" / "h4_ablation.json"))
    ap.add_argument("--study", default=str(ROOT / "docs" / "data" / "prior_corpus_study.json"))
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    budgets = [int(b) for b in args.budgets.split(",")]
    spacings = [float(s) for s in args.spacings.split(",")]
    t0 = time.time()

    footprint = footprint_from_sample(ROOT / "data" / "sample_submission.tif")
    catalogue = read_labels(ROOT / "data" / "labels.tif")
    with rasterio.open(ROOT / "data" / "external" / "sgmc" / "derived_sgmc_faults_100m_u8.tif") as s:
        sgmc = s.read(1)
    off = sgmc_off_catalogue(sgmc, catalogue, footprint)
    iso = isolated_catalogue_components(catalogue, footprint, min_separation_px=6, min_size_px=12)
    instruments = {"sgmc": dict(truth=off, mask=catalogue),
                   "iso6": dict(truth=iso["truth"], mask=iso["mask"])}
    allowed = footprint & ~catalogue

    rtp, rtp_ok = read_band(ROOT / "data" / "training_features.tif", BANDS["rtp"])
    grav, grav_ok = read_band(ROOT / "data" / "training_features.tif", BANDS["iso_grav_anom"])
    det, det_ok = read_band(ROOT / "data" / "training_features.tif", BANDS["det_elev"])
    base, base_ok = read_band(ROOT / "data" / "training_features.tif", BANDS["depth_to_base_surf"])
    ieq, ieq_ok = read_band(ROOT / "data" / "training_features.tif", BANDS["ieq_n100a15"])
    for ok in (rtp_ok, grav_ok, det_ok, base_ok, ieq_ok):
        ok &= footprint
    corrob = {"magnetic_tilt_contact": np.exp(-3.0 * np.abs(tilt_angle(rtp, rtp_ok))),
              "gravity_hgm": horizontal_gradient_magnitude(grav, grav_ok),
              "topographic_curvature": curvature_evidence(det, det_ok),
              "basement_hgm": horizontal_gradient_magnitude(base, base_ok),
              "seismicity": np.where(ieq_ok, ieq, np.nan)}
    mag = load_cloud(WORK / "cloud_rtp.npz")
    gravc = load_cloud(WORK / "cloud_iso_grav_anom.npz")

    study = json.loads(Path(args.study).read_text())
    obs = [dict(id=r["key"], n_px=float(r["emitted_px"]),
                h_sur=float(r["instruments"]["sgmc_off_catalogue"]["full_mass"]["mean_credit"]),
                kappa_sur=float(r["instruments"]["sgmc_off_catalogue"]["full_mass"]["kappa"]),
                official=float(r["score_owner_reported"]))
           for r in study["priors"] if r.get("use_for_calibration", True)]
    fit = fit_model(obs)
    print(f"calibration gamma={fit.gamma:.4f} G={fit.g_hidden:.0f} loo_rmse={fit.loo_rmse:.4f} "
          f"loo_rho={fit.loo_spearman:.3f}", flush=True)

    report = {"generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "calibration": dict(gamma=fit.gamma, g_hidden=fit.g_hidden, rmse=fit.rmse,
                                  r2=fit.r2, loo_rmse=fit.loo_rmse, loo_spearman=fit.loo_spearman),
              "budgets": budgets, "spacings": spacings,
              "instruments": {k: dict(truth_px=int(v["truth"].sum())) for k, v in instruments.items()},
              "configs": {}}
    names = [n for n in CONFIGS if not args.only or n in args.only.split(",")]
    for name in names:
        cfg = H4Config(**CONFIGS[name])
        rank, stats = assemble_ranking(mag, gravc, corrob, footprint, catalogue, cfg)
        entry = dict(overrides=CONFIGS[name], field_stats={k: v for k, v in stats.items()
                                                           if k != "config"}, lattices={})
        for sp in spacings:
            dots = greedy_min_separation(rank, allowed, n_max=200_000, min_separation_px=sp)
            yy, xx = np.nonzero(dots)
            vals = rank[yy, xx]
            order = np.argsort(-vals, kind="stable")
            yy, xx = yy[order], xx[order]
            curve = []
            for n in budgets:
                k = int(min(n, yy.size))
                m = np.zeros(rank.shape, dtype=bool)
                if k:
                    m[yy[:k], xx[:k]] = True
                row = dict(n_px=k)
                for iname, icfg in instruments.items():
                    r = binary_dti_kdtree(m, icfg["truth"], icfg["mask"], footprint)
                    row[iname] = dict(h_sur=round(r["mean_credit"], 5), kappa_sur=round(r["kappa"], 5),
                                      dti=round(r["dti"], 5))
                row["modelled_dti"] = round(predict_score(row["n_px"], row["sgmc"]["h_sur"],
                                                          row["sgmc"]["kappa_sur"], fit.gamma,
                                                          fit.g_hidden), 5)
                curve.append(row)
            best = max(curve, key=lambda r: r["modelled_dti"])
            entry["lattices"][f"{sp:.4f}"] = dict(capacity=int(yy.size), curve=curve,
                                                  best=dict(best, min_separation_px=sp))
        report["configs"][name] = entry
        b = entry["lattices"][f"{spacings[0]:.4f}"]["best"]
        print(f"{name:26s} sp={spacings[0]:.2f} n={b['n_px']:6d} sgmc_h={b['sgmc']['h_sur']:.4f} "
              f"iso6_h={b['iso6']['h_sur']:.4f} modelled={b['modelled_dti']:.4f} "
              f"capacity={entry['lattices'][f'{spacings[0]:.4f}']['capacity']} ({time.time()-t0:.0f}s)",
              flush=True)
    best_cfg = max(report["configs"].items(),
                   key=lambda kv: max(x["best"]["modelled_dti"] for x in kv[1]["lattices"].values()))
    report["winner"] = dict(name=best_cfg[0], overrides=best_cfg[1]["overrides"])
    Path(args.out).write_text(json.dumps(report, indent=1))
    print(f"winner: {best_cfg[0]}  -> {args.out} ({time.time()-t0:.0f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
