#!/usr/bin/env python3
"""Stage 1 of H4: build and cache the multi-height, multi-window Euler clouds.

Reads the official 19-band feature stack (sha256-pinned), deconvolves the two
independent potential fields (reduced-to-pole magnetics and isostatic gravity)
at three continuation heights and three window sizes with structural index 0,
and caches the resulting depth-labelled point clouds under ``work/h4/`` so the
emission stage can be re-run without repeating the inversion.

Nothing here writes a submission.  Output: ``work/h4/cloud_<layer>.npz`` plus
``work/h4/cloud_stats.json``.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe40 import BANDS  # noqa: E402
from gemsdoe40.euler_h4 import StackSpec, build_cloud  # noqa: E402
from gemsdoe40.grid import footprint_from_sample, read_band  # noqa: E402

LAYERS = {"rtp": 2, "iso_grav_anom": 13}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", default=str(ROOT / "data" / "training_features.tif"))
    ap.add_argument("--sample", default=str(ROOT / "data" / "sample_submission.tif"))
    ap.add_argument("--out", default=str(ROOT / "work" / "h4"))
    ap.add_argument("--heights", default="0,500,1500")
    ap.add_argument("--windows", default="6,8,12")
    ap.add_argument("--stride", type=int, default=4)
    ap.add_argument("--analytic-percentile", type=float, default=72.0)
    ap.add_argument("--max-rel-se", type=float, default=0.22)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    footprint = footprint_from_sample(args.sample)
    spec = StackSpec(
        layers=tuple(LAYERS), heights_m=tuple(float(v) for v in args.heights.split(",")),
        windows_px=tuple(int(v) for v in args.windows.split(",")), stride_px=args.stride,
        structural_index=0.0, analytic_percentile=args.analytic_percentile,
        max_rel_se=args.max_rel_se,
    )
    stats = {"spec": spec.__dict__, "layers": {}}
    for layer, band in LAYERS.items():
        t0 = time.time()
        field, valid = read_band(args.features, band)
        valid = valid & footprint
        print(f"[{layer}] band {band} ({BANDS.get(layer)}) valid {int(valid.sum()):,}", flush=True)
        merged, persisted = build_cloud(layer, field, valid, spec)
        print(f"[{layer}] merged {len(merged):,}  persisted {len(persisted):,}  "
              f"({time.time() - t0:.1f}s)", flush=True)
        np.savez_compressed(
            out / f"cloud_{layer}.npz",
            row=merged.row, col=merged.col, depth_m=merged.depth_m, rel_se=merged.rel_se,
            analytic=merged.analytic, cond=merged.cond,
            p_row=persisted.row, p_col=persisted.col, p_depth_m=persisted.depth_m,
            p_rel_se=persisted.rel_se, p_analytic=persisted.analytic,
        )
        stats["layers"][layer] = {
            "band": band, "merged": int(len(merged)), "persisted": int(len(persisted)),
            "merged_stats": merged.stats, "persisted_stats": persisted.stats,
            "depth_median_persisted_m": (float(np.median(persisted.depth_m)) if len(persisted) else None),
        }
    (out / "cloud_stats.json").write_text(json.dumps(stats, indent=2))
    print(f"wrote {out}/cloud_stats.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
