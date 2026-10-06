#!/usr/bin/env python3
"""H45 — Euler deconvolution depth-clustering candidate (generation only).

Pipeline
--------
1. Read the hash-pinned competition rasters (``scripts/acquire_data.py``).
2. Run Euler deconvolution for a preregistered lattice of families
   (field x structural index x window size) using ``gemsdoe40.euler_h45``.
3. Weight every depth-labelled solution by shallowness, Euler-equation misfit, conditional
   depth precision, local depth consistency and cross-scale corroboration.
4. Convert the cloud into a continuous raster by kernel-density estimation.
5. Write to the ignored ``work/`` directory.  Nothing here reads a prior submission, a
   holdout truth raster or a proxy catalogue, and nothing here decides whether the file is
   publishable: that is ``scripts/audit_h45.py``.

Usage
-----
    . .venv/bin/activate
    python scripts/run_h45_euler.py [--out work/h45] [--fast]
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import rasterio

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gemsdoe40.euler_h45 import Family, SolveConfig, run_family        # noqa: E402
from gemsdoe40.depthcluster_h45 import (ClusterConfig, cluster_weights,  # noqa: E402
                                        cross_family_weight, kde, normalise, splat)

ROOT = Path(__file__).resolve().parents[1]

# --------------------------------------------------------------------------------------
# Preregistered design.  Structural indices come from Reid & Thurston (2014) Table 1:
#   magnetic  finite contact/fault = 0 (Reid et al. 1990 eq. 2, with the free offset A)
#   magnetic  thin sheet edge      = 1
#   gravity   thin sheet edge      = 0
# and, because applying n vertical derivatives lowers the effective index by n, running the
# Euler equation with SI = 0 on the first vertical derivative of gravity is the same model
# as SI = 1 on gravity, i.e. the "thin bed fault" row of the same table.
#
# The corrected gravity value for a *finite* contact/fault is SI = -1, but Reid & Thurston
# warn that it "requires a more generalized formulation"; it is run here as a reported
# diagnostic family so that the claim can be tested rather than assumed.
# --------------------------------------------------------------------------------------
FAMILIES = [
    # gravity: the same physical field at three indices of the corrected table
    Family("iso_grav_anom", "gravity", si=0.0, window=15, stride=3),
    Family("iso_grav_anom", "gravity", si=0.0, window=25, stride=3),
    Family("iso_grav_anom", "gravity", si=0.0, window=9, stride=3),
    Family("iso_grav_anom", "gravity", si=0.0, window=15, stride=3, derivative_order=1),
    Family("iso_grav_anom", "gravity", si=-1.0, window=15, stride=3),
    # magnetics: contact (SI=0, Reid et al. eq. 2) and thin sheet edge (SI=1)
    Family("tmi", "magnetic", si=0.0, window=9, stride=3),
    Family("tmi", "magnetic", si=0.0, window=15, stride=3),
    Family("tmi", "magnetic", si=0.0, window=25, stride=3),
    Family("tmi", "magnetic", si=1.0, window=15, stride=3),
    Family("rtp", "magnetic", si=0.0, window=15, stride=3),
]

FAST_FAMILIES = [
    Family("iso_grav_anom", "gravity", si=0.0, window=15, stride=4),
    Family("iso_grav_anom", "gravity", si=0.0, window=25, stride=4),
    Family("iso_grav_anom", "gravity", si=0.0, window=15, stride=4, derivative_order=1),
    Family("iso_grav_anom", "gravity", si=-1.0, window=15, stride=4),
    Family("tmi", "magnetic", si=0.0, window=9, stride=4),
    Family("tmi", "magnetic", si=0.0, window=15, stride=4),
    Family("tmi", "magnetic", si=1.0, window=15, stride=4),
]

SOLVE = SolveConfig(
    dx=100.0,
    min_along_strike_power=0.25,
    min_depth_m=0.0,
    max_depth_m=6000.0,
    max_rel_residual=0.30,
    max_depth_se_frac=0.50,
    max_abs_offset_frac=0.60,
    min_window_valid_frac=0.95,
    edge_suppress_cells=6,
)

CLUSTER = ClusterConfig(
    depth_scale_m=900.0,
    residual_scale=0.30,
    depth_se_frac=0.35,
    neighbour_radius_m=700.0,
    scatter_tol_m=400.0,
    min_neighbours=2,
    cross_scale_bonus=0.75,
    kde_sigma_px=2.0,
    kde_truncate=4.0,
)


def load_inputs():
    with rasterio.open(ROOT / "data/sample_submission.tif") as ds:
        template = ds.read(1)
        profile = ds.profile
    footprint = np.isfinite(template)
    with rasterio.open(ROOT / "data/training_features.tif") as ds:
        names = [d.split(" - ")[0] for d in ds.descriptions]
        bands = {n: ds.read(i + 1) for i, n in enumerate(names)}
    return footprint, bands, profile


def valid_mask(bands: dict, footprint: np.ndarray) -> np.ndarray:
    mask = footprint.copy()
    for arr in bands.values():
        mask &= np.isfinite(arr) & (np.asarray(arr) > -1e30)
    return mask


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=ROOT / "work" / "h45")
    parser.add_argument("--fast", action="store_true")
    parser.add_argument("--stride", type=int, default=None)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    footprint, bands, profile = load_inputs()
    valid = valid_mask(bands, footprint)
    print(f"footprint {int(footprint.sum())}  valid {int(valid.sum())}", flush=True)

    families = FAST_FAMILIES if args.fast else FAMILIES
    if args.stride:
        families = [Family(f.name, f.kind, f.si, f.window, args.stride,
                           f.upward_m, f.derivative_order) for f in families]

    clouds = []
    summary = []
    for fam in families:
        t = time.time()
        cloud = run_family(bands[fam.name], valid, fam, SOLVE)
        clouds.append(cloud)
        z = cloud.depth
        row = dict(family=fam.name, kind=fam.kind, si=fam.si, window=fam.window,
                   stride=fam.stride, derivative_order=fam.derivative_order,
                   solutions=int(len(cloud)),
                   depth_median=float(np.median(z)) if len(z) else None,
                   depth_p10=float(np.percentile(z, 10)) if len(z) else None,
                   depth_p90=float(np.percentile(z, 90)) if len(z) else None,
                   residual_median=float(np.median(cloud.rel_residual)) if len(z) else None,
                   seconds=round(time.time() - t, 1))
        summary.append(row)
        print("  " + json.dumps(row), flush=True)

    cluster_weights(clouds, CLUSTER, dx=SOLVE.dx)
    cross_family_weight(clouds, CLUSTER, dx=SOLVE.dx)

    # cache the solution cloud so that evaluation sweeps do not re-run the deconvolution
    cache = {}
    for i, cloud in enumerate(clouds):
        tag = f"{i:02d}_{cloud.family}_si{cloud.si}_w{cloud.window}_d{cloud.meta.get('derivative_order', 0)}"
        cache[tag] = dict(col=cloud.col, row=cloud.row, depth=cloud.depth,
                          depth_se=cloud.depth_se, rel_residual=cloud.rel_residual,
                          weight=cloud.meta.get("weight", np.ones(len(cloud))),
                          family=cloud.family, kind=cloud.kind, si=cloud.si,
                          window=cloud.window)
    np.savez_compressed(args.out / "clouds.npz",
                        **{f"{k}__{f}": v for k, d in cache.items() for f, v in d.items()})

    shape = valid.shape
    per_family = {}
    for cloud in clouds:
        raw = splat(cloud, shape, dx=SOLVE.dx)
        field = kde(raw, CLUSTER.kde_sigma_px, CLUSTER.kde_truncate)
        key = f"{cloud.family}|si{cloud.si}|w{cloud.window}|d{cloud.meta.get('derivative_order', 0)}"
        per_family[key] = field
        print(f"    {key}: mass={field.sum():.1f} support={int((field > 0).sum())}", flush=True)
        np.save(args.out / f"field_{key.replace('|', '_')}.npy", field.astype(np.float32))

    total = np.zeros(shape, dtype=np.float64)
    for key, field in per_family.items():
        pos = field[valid & (field > 0)]
        scale = float(np.percentile(pos, 99.5)) if pos.size else 0.0
        if scale > 0:
            total += np.clip(field / scale, 0.0, 1.0)
    total = normalise(total, valid, mode="max")

    np.save(args.out / "combined_field.npy", total.astype(np.float32))
    np.save(args.out / "valid_mask.npy", valid)
    receipt = dict(
        generated_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        families=summary,
        solve=SOLVE.__dict__, cluster=CLUSTER.__dict__,
        combined_mass=float(np.nansum(total[valid])),
        combined_support=int((total > 0).sum()),
        elapsed_s=round(time.time() - t0, 1),
        provenance="hash-pinned public mirrors; not an authenticated DrivenData download",
    )
    (args.out / "generation_receipt.json").write_text(json.dumps(receipt, indent=1))
    print(json.dumps(receipt["families"], indent=1))
    print(f"combined mass {receipt['combined_mass']:.1f}  support {receipt['combined_support']}"
          f"  elapsed {receipt['elapsed_s']}s")


if __name__ == "__main__":
    main()
