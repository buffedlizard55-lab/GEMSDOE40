#!/usr/bin/env python3
"""Measure, for every competition layer, how strongly its high values sit on real faults.

Motivation
----------
Before designing any detector we need to know which of the 19 supplied GeoDAWN layers
carry *any* usable near-surface fault signal at the 100 m grid, and which are dominated by
deep/regional structure.  This is a purely descriptive measurement against the public
USGS/INGENIOUS label raster (``data/labels.tif``), reported two ways:

* pooled: every valid pixel in the footprint,
* blocked: a K x L spatial partition, pooled *within* each block and then averaged, so a
  single dense fault swarm cannot drive the whole result.

A layer with blocked lift ~ 1.0 carries no spatially generalising fault signal.

No hidden label, proxy model or prior submission is used.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import rasterio
from scipy.ndimage import distance_transform_edt, uniform_filter, median_filter

# ----------------------------------------------------------------------------------
# Diagnostic transforms.  ``abs`` variants catch edges irrespective of sign; the
# percentile rank is what actually gets measured, so only the ordering matters.
# ----------------------------------------------------------------------------------
def _rank(a: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Percentile rank in [0,1] over valid pixels; invalid -> 0."""
    out = np.zeros(a.shape, dtype=np.float32)
    v = a[valid]
    if v.size == 0:
        return out
    order = np.argsort(v, kind="stable")
    ranks = np.empty(v.size, dtype=np.float64)
    ranks[order] = np.arange(1, v.size + 1, dtype=np.float64)
    out[valid] = (ranks / v.size).astype(np.float32)
    return out


def _grad_mag(a: np.ndarray, valid: np.ndarray) -> np.ndarray:
    b = np.where(valid, a, 0.0)
    gy, gx = np.gradient(b)
    return np.hypot(gy, gx)


def _laplace(a: np.ndarray, valid: np.ndarray) -> np.ndarray:
    from scipy.ndimage import laplace
    return np.abs(laplace(np.where(valid, a, 0.0)))


def _ridge(a: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Local-maximum-in-the-gradient-direction thinning: keep |grad| where it is a
    transverse maximum.  Cheap, deterministic, no threshold-on-amplitude."""
    g = _grad_mag(a, valid)
    from scipy.ndimage import maximum_filter
    return np.where(g >= maximum_filter(g, size=3), g, 0.0)


def build_transforms(stack: dict[str, np.ndarray], valid: np.ndarray) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    for name, band in stack.items():
        out[name] = _rank(np.abs(np.asarray(band, dtype=np.float64)), valid)
        out[f"{name}|gradmag"] = _rank(_grad_mag(np.asarray(band, dtype=np.float64), valid), valid)
    # curvature-like combinations that several sibling repos used, recomputed here
    tmi = stack.get("tmi")
    if tmi is not None:
        tmi = np.asarray(tmi, dtype=np.float64)
        out["tmi|laplace"] = _rank(_laplace(tmi, valid), valid)
        out["tmi|ridge"] = _rank(_ridge(tmi, valid), valid)
    elev = stack.get("det_elev")
    if elev is not None:
        elev = np.asarray(elev, dtype=np.float64)
        out["det_elev|ridge"] = _rank(_ridge(elev, valid), valid)
    return out


def lift_at_top(field: np.ndarray, target: np.ndarray, valid: np.ndarray, budget: int) -> float:
    """Fraction of the top-``budget`` pixels that lie within 300 m of ``target``,
    divided by the footprint-wide base rate."""
    f = np.where(valid, field, -np.inf)
    flat = f.ravel()
    k = min(budget, int(valid.sum()))
    idx = np.argpartition(flat, -k)[-k:]
    sel = np.zeros(valid.size, dtype=bool)
    sel[idx] = True
    sel &= valid.ravel()
    base = target.sum() / max(int(valid.sum()), 1)
    hit = target.ravel()[sel].mean() if sel.any() else 0.0
    return float(hit / base) if base > 0 else float("nan")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path("work/layer_skill.json"))
    parser.add_argument("--budget", type=int, default=40000)
    parser.add_argument("--rows", type=int, default=6)
    parser.add_argument("--cols", type=int, default=5)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    with rasterio.open("data/sample_submission.tif") as ds:
        template = ds.read(1)
    footprint = np.isfinite(template)
    with rasterio.open("data/training_features.tif") as ds:
        descs = list(ds.descriptions)
        stack = {d.split(" - ")[0]: ds.read(i + 1) for i, d in enumerate(descs)}
    with rasterio.open("data/labels.tif") as ds:
        labels = ds.read(1)
    known = (labels == 1) & footprint

    # one validity mask consistent across bands (finite and not the -3.4e38 sentinel)
    valid = footprint.copy()
    for band in stack.values():
        valid &= np.isfinite(band) & (np.asarray(band) > -1e30)
    # shrink by 3 px so the 300 m target halo never runs off the grid
    from scipy.ndimage import binary_erosion
    valid &= binary_erosion(footprint, np.ones((7, 7), bool))

    d = distance_transform_edt(~known)
    target = (d <= 3.0) & valid                      # within 300 m of a mapped fault
    halo = (d <= 3.0)

    transforms = build_transforms(stack, valid)

    H, W = valid.shape
    by, bx = H // args.rows, W // args.cols
    blocks = [(r * by, min(H, (r + 1) * by), c * bx, min(W, (c + 1) * bx))
              for r in range(args.rows) for c in range(args.cols)]

    results = {}
    for name, field in transforms.items():
        pooled = lift_at_top(field, halo & valid, valid, args.budget)
        per_block = []
        for y0, y1, x0, x1 in blocks:
            sub_valid = valid[y0:y1, x0:x1]
            sub_halo = halo[y0:y1, x0:x1]
            if sub_valid.sum() < 1000 or sub_halo.sum() < 20:
                continue
            budget_b = max(20, int(args.budget * sub_valid.sum() / valid.sum()))
            per_block.append(lift_at_top(field[y0:y1, x0:x1], sub_halo, sub_valid, budget_b))
        results[name] = {
            "lift_pooled": pooled,
            "lift_blocked_mean": float(np.mean(per_block)) if per_block else float("nan"),
            "lift_blocked_median": float(np.median(per_block)) if per_block else float("nan"),
            "blocks_used": len(per_block),
        }

    order = sorted(results, key=lambda k: -(results[k]["lift_blocked_mean"] if np.isfinite(results[k]["lift_blocked_mean"]) else -9))
    print(f"{'transform':28s} {'pooled':>7} {'blocked':>8} {'median':>7} {'blocks':>6}")
    for k in order:
        r = results[k]
        print(f"{k:28s} {r['lift_pooled']:7.3f} {r['lift_blocked_mean']:8.3f} "
              f"{r['lift_blocked_median']:7.3f} {r['blocks_used']:6d}")

    args.out.write_text(json.dumps({
        "evidence_class": "descriptive layer diagnostic against the public label raster",
        "budget_px": args.budget, "block_grid": [args.rows, args.cols],
        "valid_px": int(valid.sum()), "target_px": int(target.sum()),
        "base_rate": float(halo[valid].sum() / valid.sum()),
        "results": results,
    }, indent=1))


if __name__ == "__main__":
    main()
