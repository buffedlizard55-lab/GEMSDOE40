#!/usr/bin/env python3
"""Rank prior rasters on the registered proxy instrument at *equal emitted mass*.

Registered instrument (docs/research/hypotheses.md, frozen 2026-10-05):
  truth   = off-catalogue pixels of the owner-derived SGMC fault mirror
  holdout = 4 rows x 6 columns contiguous blocks, three-cell guard on internal edges
  metric  = official distance-weighted Tversky index, alpha=0.2, beta=0.8, R=3 px
  emission = each raster thinned to an equal mass of 45,000 dots by 3 px
             non-maximum suppression with the exact catalogue plus a 2 px buffer
             excluded, so the comparison is size-neutral

PROVENANCE CAVEAT: the pinned proxy has SHA-256
``26d142c4c93282cd94f6950ab96f22aeff59fbbea523d43d662e76fa1b161b5c`` and lived at
``/tmp/proxy-data/`` which does not survive between sandbox turns.  The raster
available here, ``work/ext/derived_sgmc_faults_100m_u8.tif``, has SHA-256
``643cbe99...`` and is therefore NOT the pinned file.  Scores from this script
are a *diagnostic ranking only* and must not be quoted as the frozen promotion
gate result.
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

from gemsdoe40.emit import nms_dots  # noqa: E402
from gemsdoe40.research_metric import spatial_block_components  # noqa: E402

PROXY = ROOT / "work" / "ext" / "derived_sgmc_faults_100m_u8.tif"
PROXY_SHA256_ACTUAL = "643cbe99..."
PROXY_SHA256_PINNED = "26d142c4c93282cd94f6950ab96f22aeff59fbbea523d43d662e76fa1b161b5c"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, default=ROOT / "prior_cache")
    ap.add_argument("--mass", type=int, default=45000)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", type=Path, default=ROOT / "docs" / "data" / "proxy_ranking.json")
    args = ap.parse_args()

    with rasterio.open(ROOT / "data" / "sample_submission.tif") as src:
        foot = np.isfinite(src.read(1))
        shape, transform, crs = src.shape, src.transform, src.crs
    with rasterio.open(ROOT / "data" / "labels.tif") as src:
        cat = src.read(1) == 1
    with rasterio.open(PROXY) as src:
        proxy = src.read(1) > 0
    truth = proxy & ~cat
    del proxy

    # Equal-mass random controls (three draws) for reference at this mass.
    controls = []
    for seed in (11, 12, 13):
        rng = np.random.default_rng(seed)
        rand = np.zeros(shape, np.float32)
        rand.ravel()[rng.choice(np.flatnonzero(foot.ravel()), size=args.mass, replace=False)] = 1.0
        comp, blocks = spatial_block_components(rand, truth, foot, n_rows=4, n_cols=6, guard=3)
        controls.append({"seed": seed, "dti": comp.score, "tpw": comp.tp, "fpw": comp.fp,
                         "blocks": blocks})
    control_mean = float(np.mean([c["dti"] for c in controls]))

    files = sorted(args.cache.glob("*.tif"))
    if args.limit:
        files = files[: args.limit]
    rows = []
    t0 = time.time()
    for i, path in enumerate(files, start=1):
        try:
            with rasterio.open(path) as src:
                if src.shape != shape or src.crs != crs or src.transform != transform:
                    rows.append({"file": path.name, "skip": "grid mismatch"})
                    continue
                arr = np.where(np.isfinite(src.read(1)), src.read(1), 0.0).astype(np.float64)
        except Exception as exc:
            rows.append({"file": path.name, "skip": f"read error: {exc}"})
            continue
        dots = nms_dots(arr, foot, min_separation_px=3, budget=args.mass,
                        catalogue=cat, catalogue_buffer_px=2)
        comp, blocks = spatial_block_components(dots, truth, foot, n_rows=4, n_cols=6, guard=3)
        n = int((dots > 0).sum())
        rows.append({"file": path.name, "emitted_px": n, "dti": comp.score,
                     "tpw": comp.tp, "fpw": comp.fp, "credit_per_dot": float(comp.tp / max(n, 1)),
                     "blocks": blocks})
        if i % 20 == 0:
            print(f"{i}/{len(files)} ({time.time()-t0:.0f}s) best so far "
                  f"{max((r['dti'] for r in rows if 'dti' in r), default=0):.4f}", flush=True)

    scored = [r for r in rows if "dti" in r]
    scored.sort(key=lambda r: -r["dti"])
    payload = {
        "instrument": "off-catalogue SGMC mirror, 4x6 blocks, guard 3, official DTI "
                      "alpha=0.2 beta=0.8 R=3px, equal-mass 3 px NMS emission",
        "mass": args.mass,
        "provenance_caveat": {
            "pinned_sha256": PROXY_SHA256_PINNED,
            "actual_raster_sha256_prefix": PROXY_SHA256_ACTUAL,
            "note": "the pinned raster is unavailable in this sandbox (/tmp is not persisted); "
                    "these numbers are a diagnostic ranking, not the frozen promotion-gate result",
        },
        "random_control": {"draws": [c["dti"] for c in controls], "mean": control_mean},
        "n_scored": len(scored),
        "best": scored[0]["file"] if scored else None,
        "rows": rows,
    }
    args.out.write_text(json.dumps(payload, indent=2, default=str) + "\n")
    print(f"random control mean {control_mean:.4f} (draws {[round(c['dti'],4) for c in controls]})")
    for r in scored[:12]:
        print(f"  {r['dti']:.4f}  {r['file'][:64]}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
