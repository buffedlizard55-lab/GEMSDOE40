#!/usr/bin/env python
"""ARCHIVED ONLY — opt-in historical reproduction, not current submission advice.

Run Euler deconvolution (SI = 0, fault-like contact) on the competition potential-field layers.

Usage
-----
    python scripts/run_euler.py --layer rtp --window 11 --out work/euler_rtp_w11.npz

Layers are the organiser's own band names from ``training_features.tif``
(``src/gems40/layers.py`` records the band numbers and their official descriptions).

Depth-clustering (``cluster_weights``) turns the raw per-window solution cloud into
shallowness x tightness weights; ``kde_field`` turns the weighted cloud into the continuous
solution-density raster that the submission is built from.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gems40.euler import (cluster_weights, euler_solutions, fill_invalid, fourier_gradients,
                          kde_field)
from gems40.grid import footprint_from_sample, read_band, read_profile
from gems40.layers import BANDS


def main() -> int:
    from gemsdoe40.legacy_guard import require_legacy_opt_in
    require_legacy_opt_in()
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data", help="directory holding training_features.tif")
    ap.add_argument("--layer", default="rtp", choices=sorted(BANDS))
    ap.add_argument("--si", type=float, default=0.0, help="structural index (0 = fault/contact)")
    ap.add_argument("--window", type=int, default=11, help="Euler window, full-resolution px")
    ap.add_argument("--stride", type=int, default=2)
    ap.add_argument("--upward", type=float, default=200.0, help="upward continuation, m")
    ap.add_argument("--lowpass", type=float, default=600.0, help="low-pass cut-off wavelength, m")
    ap.add_argument("--sigma-kde", type=float, default=1.5, help="KDE bandwidth, px")
    ap.add_argument("--out", default=None)
    ap.add_argument("--diag", default=None)
    args = ap.parse_args()

    t0 = time.time()
    ddir = Path(args.data)
    feats = ddir / "training_features.tif"
    sample = ddir / "example_submission.tif"
    if not sample.exists():
        sample = ddir / "sample_submission.tif"

    prof = read_profile(feats)
    foot = footprint_from_sample(sample)
    band = BANDS[args.layer]
    raw = read_band(feats, band["band"])
    grid, valid = fill_invalid(raw, ~np.isfinite(raw))
    print(f"[layers] {args.layer}: band {band['band']} ({band['description']})")
    print(f"[grid]   shape={prof['shape']} crs={prof['crs']} transform={prof['transform']}")

    tx, ty, tz = fourier_gradients(grid, cell_m=100.0,
                                   upward_continuation_m=args.upward,
                                   lowpass_wavelength_m=args.lowpass)
    print(f"[deriv]  |Tx| p99={np.percentile(np.abs(tx), 99):.3e} "
          f"|Tz| p99={np.percentile(np.abs(tz), 99):.3e}")

    import rasterio
    with rasterio.open(feats) as ds:
        transform = ds.transform
    sol = euler_solutions(grid, tx, ty, tz, transform, si=args.si, window_px=args.window,
                          stride=args.stride, valid=valid & foot, field=args.layer)
    print(f"[euler]  {len(sol)} window solutions ({time.time() - t0:.1f}s)")
    if len(sol) == 0:
        print("ERROR: no solutions"); return 2

    w = cluster_weights(sol, grid_shape=grid.shape)
    wt = w["weight"]
    print(f"[cluster] accepted={int((wt > 0).sum())} ({100.0 * (wt > 0).mean():.2f}%) "
          f"median depth accepted={np.median(sol.depth_m[wt > 0]):.0f} m "
          f"median n_neighbours={np.median(w['n_neighbours'][wt > 0]):.0f}")

    field = kde_field(grid.shape, sol.row[wt > 0], sol.col[wt > 0], wt[wt > 0],
                      sigma_px=args.sigma_kde)
    field = np.where(foot & valid, field, np.nan)

    out = Path(args.out or f"work/euler_{args.layer}_w{args.window}_si{args.si:g}.npz")
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out,
        field=field.astype(np.float32),
        row=sol.row.astype(np.int16), col=sol.col.astype(np.int16),
        depth_m=sol.depth_m.astype(np.float32), weight=wt.astype(np.float32),
        amp=sol.amp.astype(np.float64), rms=sol.rms.astype(np.float64),
        offset_px=sol.offset_px.astype(np.float32),
        amp_gate=w["amp_gate"], keep=w["keep"],
        latched=w["amp_gate"].astype(np.int8) * 100 + w["keep"].astype(np.int8),
        layer=args.layer, si=args.si, window=args.window, stride=args.stride,
        upward=args.upward, lowpass=args.lowpass, sigma_kde=args.sigma_kde)
    diag = dict(
        layer=args.layer, band=band["band"], band_name=band["band_name"], si=args.si,
        window_px=args.window, stride=args.stride, upward_m=args.upward,
        lowpass_m=args.lowpass, sigma_kde_px=args.sigma_kde,
        n_solutions=int(len(sol)), n_accepted=int((wt > 0).sum()),
        accepted_fraction=float((wt > 0).mean()),
        depth_percentiles_accepted={str(p): float(np.percentile(sol.depth_m[wt > 0], p))
                                    for p in (5, 25, 50, 75, 95)} if (wt > 0).any() else {},
        amp_percentiles={str(p): float(np.percentile(sol.amp, p)) for p in (50, 90, 99)},
        kde_max=float(np.nanmax(field)), kde_sum=float(np.nansum(field)),
        kde_positive_px=int(np.nansum(field > 0)),
        runtime_s=round(time.time() - t0, 2),
    )
    dpath = Path(args.diag or (str(out).replace(".npz", ".json")))
    dpath.write_text(json.dumps(diag, indent=2))
    print(f"[out]    {out}  diag={dpath}")
    print(json.dumps(diag, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
