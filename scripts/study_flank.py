#!/usr/bin/env python3
"""Flank-reach sweep: how far off the catalogue should a submission sit?

The organiser masks the published USGS/INGENIOUS catalogue pixel-exactly before
scoring, and the best prior in the corpus (h33-2-b2, 0.2778) was produced by
deleting every dot within 2 px of the catalogue.  Both facts pin the emission to
a band just outside the mask, but the mask width is not published.  This script
therefore sweeps the candidate band (``lo..hi`` px from the catalogue) against a
range of assumed mask widths and reports both surrogate instruments, so the
design choice is made on measurements rather than on the two anecdotes.
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
from gemsdoe40.euler import EulerCloud  # noqa: E402
from gemsdoe40.euler_h4 import (  # noqa: E402
    coherence, depth_cluster_kde, line_response, robust_norm, tilt_angle,
)
from gemsdoe40.grid import footprint_from_sample, read_band, read_labels  # noqa: E402
from gemsdoe40.instrument import (  # noqa: E402
    binary_dti_kdtree, isolated_catalogue_components, sgmc_off_catalogue,
)
from gemsdoe40.lattice import greedy_min_separation  # noqa: E402

WORK = ROOT / "work" / "h4"


def load_cloud(path: Path) -> EulerCloud:
    z = np.load(path)
    return EulerCloud(path.stem, 0.0, -1, -1, z["p_row"], z["p_col"], z["p_depth_m"],
                      z["p_rel_se"], np.zeros_like(z["p_row"]), z["p_analytic"], {})


def band(d: np.ndarray, lo: float, hi: float, footprint: np.ndarray) -> np.ndarray:
    """Surface that decays outward from ``lo`` px and is zero beyond ``hi`` px."""
    out = np.zeros(d.shape, dtype=np.float64)
    sel = (d >= lo) & (d <= hi) & footprint
    out[sel] = 1.0 / (d[sel] - (lo - 1.0))
    return out


def score_curve(surface, allowed, instruments, valid, budgets, spacing=4.0, n_max=120_000):
    dots = greedy_min_separation(surface, allowed, n_max=n_max, min_separation_px=spacing)
    yy, xx = np.nonzero(dots)
    order = np.argsort(-surface[yy, xx], kind="stable")
    yy, xx = yy[order], xx[order]
    out = []
    for n in budgets:
        k = int(min(n, yy.size))
        m = np.zeros(surface.shape, dtype=bool)
        if k:
            m[yy[:k], xx[:k]] = True
        row = dict(n_px=k)
        for iname, cfg in instruments.items():
            r = binary_dti_kdtree(m, cfg["truth"], cfg["mask"], valid)
            row[iname] = dict(dti=round(r["dti"], 5), h_sur=round(r["mean_credit"], 5),
                              kappa_sur=round(r["kappa"], 5))
        out.append(row)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--budgets", default="10000,20000,30000,40000,60000")
    ap.add_argument("--mask-widths", default="0,1,2,3")
    ap.add_argument("--spacings", default="4.0")
    ap.add_argument("--out", default=str(ROOT / "docs" / "data" / "flank_study.json"))
    args = ap.parse_args()
    budgets = [int(b) for b in args.budgets.split(",")]
    widths = [int(w) for w in args.mask_widths.split(",")]
    spacings = [float(s) for s in args.spacings.split(",")]
    t0 = time.time()

    footprint = footprint_from_sample(ROOT / "data" / "sample_submission.tif")
    catalogue = read_labels(ROOT / "data" / "labels.tif")
    with rasterio.open(ROOT / "data" / "external" / "sgmc" / "derived_sgmc_faults_100m_u8.tif") as s:
        sgmc = s.read(1)
    off = sgmc_off_catalogue(sgmc, catalogue, footprint)
    iso6 = isolated_catalogue_components(catalogue, footprint, min_separation_px=6, min_size_px=12)
    iso12 = isolated_catalogue_components(catalogue, footprint, min_separation_px=12, min_size_px=12)
    instruments = {"sgmc": dict(truth=off, mask=catalogue),
                   "iso6": dict(truth=iso6["truth"], mask=iso6["mask"])}

    d = distance_transform_edt(~catalogue)
    rtp, rtp_ok = read_band(ROOT / "data" / "training_features.tif", BANDS["rtp"])
    mag = load_cloud(WORK / "cloud_rtp.npz")
    grav = load_cloud(WORK / "cloud_iso_grav_anom.npz")
    kde = robust_norm(depth_cluster_kde(mag, footprint.shape, sigma_px=1.6) +
                      depth_cluster_kde(grav, footprint.shape, sigma_px=1.6), footprint)
    kb = depth_cluster_kde(mag, footprint.shape, sigma_px=3.2)
    line, _ = line_response(kde, half_len=3)
    line = robust_norm(line, footprint)
    tilt = np.exp(-3.0 * np.abs(tilt_angle(rtp, rtp_ok)))
    tilt = robust_norm(np.nan_to_num(tilt, nan=0.0), footprint)
    coh = coherence(kb, sigma=2.0)
    coh = robust_norm(np.nan_to_num(coh, nan=0.0), footprint)

    bands = {"b1_6": (1.0, 6.0), "b2_6": (2.0, 6.0), "b3_6": (3.0, 6.0),
             "b2_8": (2.0, 8.0), "b3_8": (3.0, 8.0), "b2_12": (2.0, 12.0)}
    surfaces = {}
    for name, (lo, hi) in bands.items():
        f = band(d, lo, hi, footprint)
        fn = robust_norm(f, footprint)
        surfaces[name] = f
        surfaces[f"euler_x_{name}"] = np.where(footprint, kde * (0.25 + 0.75 * fn), 0.0)
        surfaces[f"eulerline_x_{name}"] = np.where(footprint, line * (0.25 + 0.75 * fn), 0.0)
        surfaces[f"corrob_x_{name}"] = np.where(footprint, np.maximum(tilt, coh) *
                                                (0.25 + 0.75 * fn), 0.0)

    report = {"generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "budgets": budgets, "mask_widths": widths, "spacings": spacings,
              "truth_px": {"sgmc": int(off.sum()), "iso6": int(iso6["truth"].sum()),
                           "iso12": int(iso12["truth"].sum())},
              "results": {}}
    for w in widths:
        allowed = footprint & ~catalogue if w == 0 else footprint & ~binary_dilation(
            catalogue, iterations=w, structure=np.ones((3, 3), bool))
        report["results"][f"mask{w}"] = {}
        for name, surf in surfaces.items():
            report["results"][f"mask{w}"][name] = {}
            for sp in spacings:
                curve = score_curve(surf, allowed, instruments, footprint, budgets, spacing=sp)
                report["results"][f"mask{w}"][name][f"{sp:.2f}"] = curve
            c = curve
            def at(n):
                m = [r for r in c if r["n_px"] == n]
                return m[0] if m else c[-1]
            a = at(min(30000, max(r["n_px"] for r in c)))
            print(f"mask{w} {name:20s} n={a['n_px']:6d} sgmc={a['sgmc']['h_sur']:.4f} "
                  f"iso6={a['iso6']['h_sur']:.4f} "
                  f"({time.time()-t0:.0f}s)", flush=True)
            report["results"][f"mask{w}"][name]["summary"] = dict(
                n=a["n_px"], sgmc_h=a["sgmc"]["h_sur"], iso6_h=a["iso6"]["h_sur"])
    out = Path(args.out)
    out.write_text(json.dumps(report, indent=1))
    print(f"wrote {out} ({time.time()-t0:.0f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
