#!/usr/bin/env python3
"""Which evidence types actually place dots near unmapped faults?

Builds every candidate evidence surface on the competition grid - the Euler
depth-cluster KDE, its line-response/coherence sharpening, the catalogue flank
ring, catalogue tip protraction, topographic curvature, magnetic tilt-angle and
gravity horizontal-gradient lineaments, the two best priors, and a uniform
random control - emits each one as a minimum-separation dot lattice at the same
budgets, and scores them on both surrogate instruments.

Output (``docs/data/evidence_study.json``): the per-surface emit curves plus the
calibrated model score, so the design decision is made on measurements.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import rasterio
from scipy.ndimage import binary_dilation, distance_transform_edt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe40 import BANDS  # noqa: E402
from gemsdoe40.calibration import fit_model, predict_score  # noqa: E402
from gemsdoe40.corpus import PRIORS, load_field, local_path  # noqa: E402
from gemsdoe40.euler import EulerCloud  # noqa: E402
from gemsdoe40.euler_h4_line import (  # noqa: E402
    H4Config, assemble_ranking, curvature_evidence, depth_cluster_kde,
    horizontal_gradient_magnitude, line_response, robust_norm, tilt_angle,
)
from gemsdoe40.grid import footprint_from_sample, read_band, read_labels  # noqa: E402
from gemsdoe40.instrument import (  # noqa: E402
    binary_dti_kdtree, isolated_catalogue_components, random_dots, sgmc_off_catalogue,
)
from gemsdoe40.lattice import greedy_min_separation  # noqa: E402

WORK = ROOT / "work" / "h4"


def load_cloud(path: Path, prefix: str = "p_") -> EulerCloud:
    z = np.load(path)
    row = z[f"{prefix}row"]
    return EulerCloud(path.stem.replace("cloud_", ""), 0.0, -1, -1, row, z[f"{prefix}col"],
                      z[f"{prefix}depth_m"], z[f"{prefix}rel_se"], np.zeros_like(row),
                      z[f"{prefix}analytic"], {})


def tip_protraction(catalogue: np.ndarray, footprint: np.ndarray, *, reach_px: int = 25,
                    orient_radius: float = 8.0) -> np.ndarray:
    """Extend every catalogue endpoint along the local strike of the mapped trace.

    Endpoints are catalogue pixels whose local catalogue density is low; the
    local strike is the leading eigenvector of the pixel covariance within
    ``orient_radius``.  The returned surface decays with extension distance, so
    the lattice visits the most likely tips first.  This encodes the structural
    observation that expert revisions of a fault map are dominated by
    along-strike continuations of mapped traces.
    """
    from scipy.spatial import cKDTree

    cat = np.asarray(catalogue, bool) & np.asarray(footprint, bool)
    yy, xx = np.nonzero(cat)
    pts = np.column_stack([yy, xx]).astype(np.float64)
    tree = cKDTree(pts)
    dens = np.array([len(tree.query_ball_point(p, 4.0)) for p in pts])
    ends = np.flatnonzero(dens <= 4)
    out = np.zeros(catalogue.shape, dtype=np.float64)
    h, w = out.shape
    for i in ends:
        nb = tree.query_ball_point(pts[i], orient_radius)
        if len(nb) < 4:
            continue
        loc = pts[np.asarray(nb, dtype=int)] - pts[i]
        cov = loc.T @ loc
        val, vec = np.linalg.eigh(cov)
        d = vec[:, -1]  # leading eigenvector = local strike
        centroid = loc.mean(axis=0)
        if np.dot(d, centroid) < 0:
            d = -d  # point away from the trace body
        for step in range(1, reach_px + 1):
            r = int(round(pts[i][0] + d[0] * step))
            c = int(round(pts[i][1] + d[1] * step))
            if 0 <= r < h and 0 <= c < w and footprint[r, c]:
                v = 1.0 / (1.0 + 0.35 * step)
                if v > out[r, c]:
                    out[r, c] = v
    return out


def flank_surface(catalogue: np.ndarray, footprint: np.ndarray, *, reach_px: int = 6) -> np.ndarray:
    """Decaying surface around the catalogue: nearest ring first."""
    d = distance_transform_edt(~np.asarray(catalogue, bool))
    v = np.where((d > 0.5) & (d <= reach_px) & footprint, 1.0 / np.maximum(d, 1.0), 0.0)
    return v


def emit_curve(surface: np.ndarray, allowed: np.ndarray, instruments: dict, valid: np.ndarray,
               budgets, spacing: float = 2.8284, n_max: int = 120_000) -> list[dict]:
    dots = greedy_min_separation(surface, allowed, n_max=n_max, min_separation_px=spacing)
    yy, xx = np.nonzero(dots)
    vals = surface[yy, xx]
    order = np.argsort(-vals, kind="stable")
    out = []
    for n in budgets:
        k = int(min(n, yy.size))
        m = np.zeros(surface.shape, dtype=bool)
        if k:
            m[yy[order[:k]], xx[order[:k]]] = True
        row = dict(n_px=k)
        for iname, cfg in instruments.items():
            r = binary_dti_kdtree(m, cfg["truth"], cfg["mask"], valid)
            row[iname] = dict(dti=r["dti"], h_sur=r["mean_credit"], kappa_sur=r["kappa"], tp=r["tp"])
        out.append(row)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--budgets", default="10000,20000,30000,40000,60000,90000")
    ap.add_argument("--spacings", default="2.8284")
    ap.add_argument("--out", default=str(ROOT / "docs" / "data" / "evidence_study.json"))
    ap.add_argument("--study", default=str(ROOT / "docs" / "data" / "prior_corpus_study.json"))
    args = ap.parse_args()
    budgets = [int(b) for b in args.budgets.split(",")]
    spacings = [float(s) for s in args.spacings.split(",")]
    t0 = time.time()

    sample = ROOT / "data" / "sample_submission.tif"
    footprint = footprint_from_sample(sample)
    catalogue = read_labels(ROOT / "data" / "labels.tif")
    with rasterio.open(ROOT / "data" / "external" / "sgmc" / "derived_sgmc_faults_100m_u8.tif") as s:
        sgmc = s.read(1)
    off = sgmc_off_catalogue(sgmc, catalogue, footprint)
    iso = isolated_catalogue_components(catalogue, footprint, min_separation_px=6, min_size_px=12)
    iso12 = isolated_catalogue_components(catalogue, footprint, min_separation_px=12,
                                          min_size_px=12)
    print(f"instrument truth px: sgmc={int(off.sum())} iso6={int(iso['truth'].sum())} "
          f"iso12={int(iso12['truth'].sum())}", flush=True)
    instruments = {"sgmc_off_catalogue": dict(truth=off, mask=catalogue),
                   "isolated_components": dict(truth=iso["truth"], mask=iso["mask"]),
                   "isolated_sep12": dict(truth=iso12["truth"], mask=iso12["mask"])}
    allowed = footprint & ~catalogue

    rtp, rtp_ok = read_band(ROOT / "data" / "training_features.tif", BANDS["rtp"])
    grav, grav_ok = read_band(ROOT / "data" / "training_features.tif", BANDS["iso_grav_anom"])
    det, det_ok = read_band(ROOT / "data" / "training_features.tif", BANDS["det_elev"])
    rtp_ok &= footprint
    grav_ok &= footprint
    det_ok &= footprint

    mag = load_cloud(WORK / "cloud_rtp.npz")
    gravc = load_cloud(WORK / "cloud_iso_grav_anom.npz")
    shape = footprint.shape

    kde_m = depth_cluster_kde(mag, shape, sigma_px=1.6)
    kde_g = depth_cluster_kde(gravc, shape, sigma_px=1.6)
    primary = kde_m + kde_g
    line, _ = line_response(robust_norm(primary, footprint), half_len=3)

    flank = flank_surface(catalogue, footprint)
    tip = tip_protraction(catalogue, footprint)
    prim_n = robust_norm(primary, footprint)
    line_n = robust_norm(line, footprint)
    corr_max = np.maximum(np.maximum(robust_norm(np.exp(-3.0 * np.abs(tilt_angle(rtp, rtp_ok))), footprint),
                                    robust_norm(np.nan_to_num(horizontal_gradient_magnitude(grav, grav_ok), nan=0.0), footprint)),
                          robust_norm(np.nan_to_num(curvature_evidence(det, det_ok), nan=0.0), footprint))
    flank_n = robust_norm(flank, footprint)
    tip_n = robust_norm(tip, footprint)
    prox_n = np.maximum(flank_n, tip_n)

    surfaces: dict[str, np.ndarray] = {
        "euler_kde_mag": np.where(footprint, kde_m, 0.0),
        "euler_kde_both": np.where(footprint, primary, 0.0),
        "euler_line_response": np.where(footprint, line, 0.0),
        "catalogue_flank_ring": flank,
        "catalogue_tip_protraction": tip,
        "topographic_curvature": np.nan_to_num(curvature_evidence(det, det_ok), nan=0.0) * footprint,
        "magnetic_tilt_contact": np.exp(-3.0 * np.abs(tilt_angle(rtp, rtp_ok))) * footprint,
        "gravity_hgm": np.nan_to_num(horizontal_gradient_magnitude(grav, grav_ok), nan=0.0) * footprint,
        "proximity_only": np.where(footprint, prox_n, 0.0),
        # hybrid: Euler depth-cluster field modulated by structural proximity
        "euler_x_prox": np.where(footprint, prim_n * (0.25 + 0.75 * prox_n), 0.0),
        "euler_x_prox_x_line": np.where(footprint, prim_n * (0.25 + 0.75 * prox_n) * (0.5 + 0.5 * line_n), 0.0),
        "euler_x_prox_x_corrob": np.where(footprint, prim_n * (0.25 + 0.75 * prox_n) * (0.5 + 0.5 * corr_max), 0.0),
        "euler_line_x_prox": np.where(footprint, line_n * (0.25 + 0.75 * prox_n), 0.0),
        "prox_x_corrob": np.where(footprint, prox_n * (0.4 + 0.6 * corr_max), 0.0),
    }
    # priors as references, re-emitted on the same lattice
    for p in PRIORS:
        if p.key in ("h27-4-r1-solo-d2-8", "h32-1-prethin-tip-euler-d2-8", "dotted-h19-5-d2-8",
                     "h16-1-topo-geophys-ridges"):
            path = local_path(p)
            if path.exists():
                surfaces[f"prior_{p.key}"] = load_field(path, footprint)

    report = {"generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "budgets": budgets, "spacings": spacings, "surfaces": {}, "controls": {}}
    for name, surf in surfaces.items():
        report["surfaces"][name] = {}
        for sp in spacings:
            curve = emit_curve(surf, allowed, instruments, footprint, budgets, spacing=sp)
            report["surfaces"][name][f"{sp:.4f}"] = curve
        c = report["surfaces"][name][f"{spacings[0]:.4f}"]
        at40 = [r for r in c if r["n_px"] == min(40000, max(x["n_px"] for x in c))]
        at40 = at40[0] if at40 else c[-1]
        best40 = at40["sgmc_off_catalogue"]["h_sur"]
        iso40 = at40["isolated_components"]["h_sur"]
        iso12 = at40["isolated_sep12"]["h_sur"]
        print(f"{name:30s} sgmc_h={best40:.4f}  iso6_h={iso40:.4f}  iso12_h={iso12:.4f} "
              f"({time.time()-t0:.0f}s)",
              flush=True)

    for n in (20000, 40000):
        hs, ks = [], []
        for seed in range(3):
            dots = random_dots(footprint, catalogue, n, seed=seed)
            for iname, cfg in instruments.items():
                r = binary_dti_kdtree(dots, cfg["truth"], cfg["mask"], footprint)
                (hs if iname == "sgmc_off_catalogue" else ks).append(r["mean_credit"])
        report["controls"][f"random_{n}"] = dict(h_sur_sgmc=float(np.mean(hs)),
                                                 h_sur_iso=float(np.mean(ks)))

    # calibrate on the study and score every surface at its best budget
    study = json.loads(Path(args.study).read_text())
    obs = [dict(id=r["key"], n_px=float(r["emitted_px"]),
                h_sur=float(r["instruments"]["sgmc_off_catalogue"]["full_mass"]["mean_credit"]),
                kappa_sur=float(r["instruments"]["sgmc_off_catalogue"]["full_mass"]["kappa"]),
                official=float(r["score_owner_reported"]))
           for r in study["priors"] if r.get("use_for_calibration", True)]
    fit = fit_model(obs)
    report["calibration"] = dict(gamma=fit.gamma, g_hidden=fit.g_hidden, rmse=fit.rmse,
                                 r2=fit.r2, loo_rmse=fit.loo_rmse, loo_spearman=fit.loo_spearman)
    for name, per_sp in report["surfaces"].items():
        for sp, curve in per_sp.items():
            preds = [predict_score(c["n_px"], c["sgmc_off_catalogue"]["h_sur"],
                                   c["sgmc_off_catalogue"]["kappa_sur"], fit.gamma, fit.g_hidden)
                     for c in curve]
            best_i = int(np.argmax(preds))
            per_sp[sp] = {"curve": curve, "best": dict(curve[best_i], modelled_dti=float(preds[best_i]))}
    out = Path(args.out)
    out.write_text(json.dumps(report, indent=1))
    print(f"wrote {out} ({time.time()-t0:.0f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
