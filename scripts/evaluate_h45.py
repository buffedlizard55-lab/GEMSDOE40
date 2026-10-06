#!/usr/bin/env python3
"""Measure H45 family skill on the non-circular blocked holdout and fit the emission.

Reads the cached solution cloud written by ``scripts/run_h45_euler.py`` and, for every
family and for preregistered combinations of them, reports

* ``dti_blocked_mean``   - the blocked-holdout distance-weighted Tversky index
* ``lift_blocked_mean``  - top-budget hit rate against held-out mapped faults / base rate
* ``catalogue_gradient`` - where the mass sits relative to the mapped catalogue

No hidden label and no prior submission is used.  See ``src/gemsdoe40/holdout_h45.py``.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe40.depthcluster_h45 import kde, splat                    # noqa: E402
from gemsdoe40.holdout_h45 import evaluate as holdout_evaluate       # noqa: E402

BUDGET = 40000


def load_clouds(path: Path) -> dict:
    z = np.load(path, allow_pickle=False)
    groups: dict[str, dict] = {}
    for key in z.files:
        tag, field_name = key.split("__", 1)
        groups.setdefault(tag, {})[field_name] = z[key]
    return groups


def group_field(groups: dict, keys: list[str], shape, sigma: float,
                depth_scale: float) -> np.ndarray:
    total = np.zeros(shape, dtype=np.float64)
    for k in keys:
        g = groups[k]
        w = g["weight"].astype(np.float64)
        if depth_scale:
            w = w * np.exp(-np.maximum(g["depth"], 0.0) / depth_scale)
        from gemsdoe40.euler_h45 import SolutionCloud
        cloud = SolutionCloud(family=str(g["family"]), kind=str(g["kind"]),
                              si=float(g["si"]), window=int(g["window"]),
                              col=g["col"], row=g["row"], depth=g["depth"],
                              depth_se=g["depth_se"], rel_residual=g["rel_residual"],
                              offset=np.zeros(len(g["col"])),
                              along_strike=np.zeros(len(g["col"])),
                              cross_strike=np.zeros(len(g["col"])),
                              strike_ratio=np.zeros(len(g["col"])),
                              n_obs=np.zeros(len(g["col"])))
        cloud.meta["weight"] = w
        total += splat(cloud, shape, dx=100.0)
    return kde(total, sigma, 4.0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, default=ROOT / "work" / "h45" / "clouds.npz")
    parser.add_argument("--out", type=Path, default=ROOT / "work" / "h45" / "family_skill.json")
    parser.add_argument("--sigma", type=float, default=2.0)
    parser.add_argument("--rows", type=int, default=6)
    parser.add_argument("--cols", type=int, default=5)
    args = parser.parse_args()

    with rasterio.open(ROOT / "data/sample_submission.tif") as ds:
        foot = np.isfinite(ds.read(1))
    with rasterio.open(ROOT / "data/labels.tif") as ds:
        lab = ds.read(1)
    valid = np.load(ROOT / "work" / "h45" / "valid_mask.npy")
    truth = (lab == 1) & foot & valid
    print(f"valid {int(valid.sum())} truth {int(truth.sum())}", flush=True)

    groups = load_clouds(args.cache)
    keys = sorted(groups)
    print("families:")
    for k in keys:
        g = groups[k]
        print(f"  {k:44s} n={len(g['col']):6d} depth_med={np.median(g['depth']) if len(g['depth']) else 0:8.0f}")

    results = {}
    for k in keys:
        if len(groups[k]["col"]) == 0:
            results[k] = dict(skipped="no solutions")
            continue
        field = group_field(groups, [k], valid.shape, args.sigma, depth_scale=0.0)
        field = np.where(valid, field, 0.0)
        m = field.max()
        if m <= 0:
            results[k] = dict(skipped="zero field")
            continue
        # evaluate at a fixed mass budget so families are compared like for like
        scaled = field / m
        scaled = scaled * (BUDGET / max(scaled.sum(), 1e-9))
        res = holdout_evaluate(scaled, truth, valid, rows=args.rows, cols=args.cols,
                               budget_px=BUDGET)
        results[k] = {kk: vv for kk, vv in res.items() if kk != "per_block"}
        print(f"  {k:44s} dti={res['dti_blocked_mean']:.5f} lift={res['lift_blocked_mean']:.3f} "
              f"sup={res['support']}", flush=True)

    # shallow-weighted and combined variants
    variants = {
        "ALL_unweighted": (keys, 0.0),
        "ALL_shallow900": (keys, 900.0),
        "ALL_shallow400": (keys, 400.0),
        "MAG_only": ([k for k in keys if "tmi" in k or "rtp" in k], 900.0),
        "GRAV_only": ([k for k in keys if "grav" in k and "si-1.0" not in k], 900.0),
    }
    for name, (ks, ds_) in variants.items():
        ks = [k for k in ks if len(groups[k]["col"]) > 0]
        if not ks:
            continue
        field = group_field(groups, ks, valid.shape, args.sigma, ds_)
        field = np.where(valid, field, 0.0)
        m = field.max()
        if m <= 0:
            continue
        scaled = field / m * (BUDGET / max((field / m).sum(), 1e-9))
        res = holdout_evaluate(scaled, truth, valid, rows=args.rows, cols=args.cols,
                               budget_px=BUDGET)
        results[name] = {kk: vv for kk, vv in res.items() if kk != "per_block"}
        print(f"  {name:44s} dti={res['dti_blocked_mean']:.5f} lift={res['lift_blocked_mean']:.3f}",
              flush=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(results, indent=1, default=float))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
