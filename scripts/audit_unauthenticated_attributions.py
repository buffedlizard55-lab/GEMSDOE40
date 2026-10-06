#!/usr/bin/env python3
"""Audit record of a PROHIBITED method: correlating emitted-dot properties with
unauthenticated, group-reported leaderboard scores.

WHY THIS FILE EXISTS.  A leaderboard row establishes an account-level public
score only; it does not identify a raster unless the organizers link the
submission file (see ``docs/data/feed-20261005.json``, which specifically
forbids attributing the 0.2778 row to ``gemsdoe32-h33-h33-2-b2-…``).  Every
pairing used here is unauthenticated.  The correlations this script computes
are therefore NOT evidence for any design decision, emission weight, threshold
or instrument judgement, and the binding rule in the register
(``docs/research/hypotheses.md``, amendment of 2026-10-06) forbids using them.

The script is retained so the trap is documented and not re-entered: it
reproduces the apparent correlations (n = 10, two-cluster structure) that were
briefly taken as a design basis and then withdrawn.  Output carries
``status = WITHDRAWN_AS_EVIDENCE``.

Do not import from this module and do not derive weights from its output.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import rasterio
from scipy import ndimage, stats
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

GROUP_REPORTED_UNAUTHENTICATED_SCORES = {
    "gemsdoe32-h33-h33-2-b2-20261004T220000Z-e5eb6e7e-nan.tif": 0.2778,
    "gems28-h27-4-r1-solo-d2-8-20261003-8acb75e1f2cc-nan.tif": 0.2708,
    "gems25-dotted-h19-5-d2-8-20261002-e56ea318af89-zeros.tif": 0.2600,
    "gems24-h25-1-dotted-h19-5-d1-5-20261002-989f59505db1-nan.tif": 0.2477,
    "gems27-topo-gap-closure-t-v2-on-d1-5-20261002-5512495c6bd1-nan.tif": 0.2449,
    "gems24-reference-h19-5-20261002-80d47e1ab2ee-nan.tif": 0.1922,
    "gems10-h16-continuation-20260927T065521077735Z-3431b83c7c.tif": 0.0461,
    "gems10-h20-dem10-scarp-thin-20260927T155223039488Z-ffc91a1686.tif": 0.0921,
    "gems10-h25-ctx-ridge-20260927T232947704150Z-6452ae1d00.tif": 0.1280,
    "gems10-h28-dotted-ridge-20260928T020256236880Z-6452ae1d00.tif": 0.1839,
}


def main() -> int:
    inventory = json.loads((ROOT / "docs" / "data" / "prior_raster_inventory.json").read_text())
    name2blob: dict[str, str] = {}
    for entry in inventory["unique_rasters"]:
        for artifact in entry["artifacts"]:
            name2blob.setdefault(Path(artifact["path"]).name, entry["git_blob_sha"])

    with rasterio.open(ROOT / "data" / "sample_submission.tif") as src:
        foot = np.isfinite(src.read(1))
    cat = rasterio.open(ROOT / "data" / "labels.tif").read(1) == 1
    dist_cat = ndimage.distance_transform_edt(~cat)
    cat_buf = dist_cat <= 3.0

    sgmc_path = ROOT / "work" / "ext" / "derived_sgmc_faults_100m_u8.tif"
    sgmc_off = None
    if sgmc_path.exists():
        with rasterio.open(sgmc_path) as src:
            sgmc_off = (src.read(1) > 0) & ~cat

    from gemsdoe40.depthcluster import hessian_lineament_response, robust_normalize, terrain_concurrence

    slope = rasterio.open(ROOT / "data" / "training_features.tif").read(19).astype(np.float64)
    with rasterio.open(ROOT / "work" / "ext" / "lidar_scarp_features_u8.tif") as src:
        upface = src.read(7).astype(np.float64)
        lidar_cov = src.read(12).astype(np.float64) > 0

    bands: dict[str, np.ndarray] = {}
    with rasterio.open(ROOT / "data" / "training_features.tif") as src:
        for b in (3, 4, 12, 16, 18, 19):
            arr = src.read(b).astype(np.float64)
            arr[np.abs(arr) > 1e29] = np.nan
            bands[src.tags(b).get("band_name", str(b))] = arr

    rows: dict[str, dict] = {}
    for name, live in GROUP_REPORTED_UNAUTHENTICATED_SCORES.items():
        blob = name2blob.get(name)
        if not blob:
            print(f"skip (no blob): {name}")
            continue
        path = ROOT / "prior_cache" / f"{blob}.tif"
        if not path.exists():
            print(f"skip (not cached): {name}")
            continue
        with rasterio.open(path) as src:
            arr = np.where(np.isfinite(src.read(1)), src.read(1), 0.0).astype(np.float64)
        pos = arr > 0
        yy, xx = np.nonzero(pos)
        n = int(pos.sum())
        if n == 0:
            continue
        props: dict[str, float] = {"live": float(live), "n_px": float(n)}
        props["median_dist_cat"] = float(np.median(dist_cat[yy, xx]))
        props["frac_on_catbuf"] = float(cat_buf[yy, xx].mean())
        if sgmc_off is not None:
            props["frac_on_sgmc_off"] = float(sgmc_off[yy, xx].mean())
        props["frac_lidar"] = float(lidar_cov[yy, xx].mean())
        props["mean_upface"] = float(np.nanmean(upface[yy, xx]))
        for key, band in bands.items():
            props[f"mean_{key}"] = float(np.nanmean(band[yy, xx]))
        sub = np.column_stack([yy, xx])
        if n > 8000:
            sub = sub[np.random.default_rng(0).choice(n, size=8000, replace=False)]
        tree = cKDTree(sub)
        dd, _ = tree.query(sub, k=2)
        props["median_nn_dots"] = float(np.median(dd[:, 1]))
        blk = (yy // 933) * 6 + (xx // 549)
        p = np.bincount(blk, minlength=24) / n
        props["block_entropy"] = float(-(p[p > 0] * np.log(p[p > 0])).sum())
        rows[name] = props

    keys = sorted({k for props in rows.values() for k in props} - {"live"})
    table = []
    for key in keys:
        xs = [rows[n].get(key, np.nan) for n in rows]
        ys = [rows[n]["live"] for n in rows]
        ok = np.isfinite(xs)
        if ok.sum() < 4:
            continue
        rho, p = stats.spearmanr(np.asarray(xs)[ok], np.asarray(ys)[ok])
        table.append({"property": key, "spearman_rho": float(rho), "p_value": float(p),
                      "n": int(ok.sum())})
    table.sort(key=lambda r: -abs(r["spearman_rho"]))

    payload = {
        "status": "WITHDRAWN_AS_EVIDENCE",
        "purpose": "record of a prohibited method: correlating emitted-dot properties with "
                   "unauthenticated group-reported leaderboard scores. Not evidence for any "
                   "design decision; see docs/research/hypotheses.md amendment 2026-10-06.",
        "prohibited": True,
        "n_files": len(rows),
        "caveats": [
            "n = 10 files from one lineage; five are variants of the same dotted-ridge family, "
            "so the correlations are partly a two-cluster contrast rather than a replicated law.",
            "Correlations between the measured properties are high (distance-to-catalogue and "
            "detrended elevation both express 'away from the mapped basin-margin faults').",
            "Scores are group-reported in the family manifests and are not re-verified by an "
            "upload from this sandbox.",
        ],
        "correlations": table,
        "per_file": rows,
    }
    out = ROOT / "docs" / "data" / "unauthenticated_attribution_audit.json"
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n")
    print("WITHDRAWN AS EVIDENCE — unauthenticated attributions; reproduce-only record.")
    for row in table:
        print(f"{row['property']:26s} rho={row['spearman_rho']:+.3f}  p={row['p_value']:.4f}  n={row['n']}")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
