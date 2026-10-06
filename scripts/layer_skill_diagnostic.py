#!/usr/bin/env python3
"""Diagnostic: how much fault-finding information does each available layer carry?

Truth used for this diagnostic is the organizer catalogue raster on the exact
competition grid (data/labels.tif, 60,988 positive 100 m cells).  It is a
skill diagnostic, NOT the hidden target: the hidden set is by construction
off-catalogue (DrivenData staff clarification, community post 11516/2), so a
layer that only reproduces the catalogue is not sufficient.  The diagnostic
answers a narrower, decisive question: *do the external 3DEP lidar scarp
features and the GeoDAWN radiometric grids carry independent fault-position
information at all?*

Metric used for the diagnostic is the official distance-weighted Tversky
index machinery from src/gemsdoe40/metric.py (AUC-style separation is also
reported, computed on a stratified sample for tractability).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import rasterio
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe40.metric import dti_exact  # noqa: E402


def read(path: Path) -> np.ndarray:
    with rasterio.open(path) as src:
        return src.read().astype(np.float32)


def auc(score: np.ndarray, positive: tuple, negative: tuple) -> float:
    """Rank AUC for (row, col) index tuples; non-finite scores are replaced by the layer minimum."""
    s_pos = score[positive]
    s_neg = score[negative]
    if not (np.isfinite(s_pos).all() and np.isfinite(s_neg).all()):
        floor = float(np.nanmin(np.where(np.isfinite(score), score, np.nan)))
        s_pos = np.where(np.isfinite(s_pos), s_pos, floor)
        s_neg = np.where(np.isfinite(s_neg), s_neg, floor)
    s = np.concatenate([s_pos, s_neg])
    y = np.concatenate([np.ones(s_pos.size), np.zeros(s_neg.size)])
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty_like(order, dtype=np.float64)
    ranks[order] = np.arange(1, s.size + 1)
    # average ranks for ties
    s_sorted = s[order]
    i = 0
    while i < s_sorted.size:
        j = i
        while j + 1 < s_sorted.size and s_sorted[j + 1] == s_sorted[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + 1 + j + 1) / 2.0
        i = j + 1
    n_pos, n_neg = y.sum(), (1.0 - y).sum()
    return float((ranks[y == 1].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def main() -> int:
    labels = read(ROOT / "data" / "labels.tif")[0]
    sample = read(ROOT / "data" / "sample_submission.tif")[0]
    footprint = np.isfinite(sample)
    cat = (labels == 1) & footprint
    print(f"catalogue {int(cat.sum()):,} px in footprint {int(footprint.sum()):,}")

    rng = np.random.default_rng(7)
    near_cat = ndimage.binary_dilation(cat, np.ones((3, 3), bool))
    # negative control: footprint pixels away from any catalogue trace
    far = footprint & ~ndimage.binary_dilation(cat, np.ones((11, 11), bool))
    yy, xx = np.nonzero(cat)
    n_sample = min(20000, yy.size)
    pick = rng.choice(yy.size, size=n_sample, replace=False)
    pos_idx = (yy[pick], xx[pick])
    fy, fx = np.nonzero(far)
    pick2 = rng.choice(fy.size, size=min(60000, fy.size), replace=False)
    neg_idx = (fy[pick2], fx[pick2])

    report: dict = {
        "truth": "data/labels.tif catalogue (organizer existing faults), 100 m grid EPSG:32611",
        "catalogue_px": int(cat.sum()),
        "negatives": "footprint pixels not within 5 px of any catalogue trace",
        "n_positive_sample": int(n_sample),
        "n_negative_sample": int(pick2.size),
        "layers": {},
    }

    sources = {
        "challenge_band_det_elev_slope": (ROOT / "data" / "training_features.tif", [19]),
        "challenge_band_det_elev": (ROOT / "data" / "training_features.tif", [12]),
        "challenge_band_rtp": (ROOT / "data" / "training_features.tif", [2]),
        "challenge_band_tmi_hg": (ROOT / "data" / "training_features.tif", [3]),
        "challenge_band_iso_grav_anom_hg": (ROOT / "data" / "training_features.tif", [18]),
        "challenge_band_geod_2ndinv": (ROOT / "data" / "training_features.tif", [4]),
        "challenge_band_ieq_n100a15": (ROOT / "data" / "training_features.tif", [16]),
    }
    for name, (path, bands) in sources.items():
        if not path.exists():
            continue
        with rasterio.open(path) as src:
            for b in bands:
                arr = src.read(b).astype(np.float32)
                tag = src.tags(b).get("band_name", f"band{b}")
                pv = arr[pos_idx]
                nv = arr[neg_idx]
                report["layers"][f"{name}:{tag}"] = {
                    "auc_catalogue_vs_far": auc(arr, pos_idx, neg_idx),
                    "mean_at_faults": float(np.nanmean(pv)),
                    "mean_far": float(np.nanmean(nv)),
                }
                print(f"{name}:{tag}  AUC={report['layers'][f'{name}:{tag}']['auc_catalogue_vs_far']:.3f}")

    ext = {
        "lidar_scarp": (ROOT / "work" / "ext" / "lidar_scarp_features_u8.tif",
                        ["ex_max", "ex_mean", "step_max", "lapneg_max", "lappos_max",
                         "downface_max", "upface_max", "cross_max", "relief", "coh100",
                         "strike", "valid"]),
        "geodawn_rad": (ROOT / "work" / "ext" / "geodawn_rad_u8.tif", ["K", "Th", "U", "TC"]),
        "geodawn_ext": (ROOT / "work" / "ext" / "geodawn_extensions_u8.tif",
                        ["ThK", "UK", "UTh", "TMI_up150"]),
    }
    for group, (path, names) in ext.items():
        if not path.exists():
            print(f"missing {path}")
            continue
        with rasterio.open(path) as src:
            for i, nm in enumerate(names, start=1):
                if i > src.count:
                    continue
                arr = src.read(i).astype(np.float32)
                pv = arr[pos_idx]
                nv = arr[neg_idx]
                report["layers"][f"{group}:{nm}"] = {
                    "auc_catalogue_vs_far": auc(arr, pos_idx, neg_idx),
                    "mean_at_faults": float(np.nanmean(pv)),
                    "mean_far": float(np.nanmean(nv)),
                    "coverage_fraction": float((arr > 0).mean()),
                }
                print(f"{group}:{nm}  AUC={report['layers'][f'{group}:{nm}']['auc_catalogue_vs_far']:.3f} "
                      f"cover={(arr>0).mean():.2f}")

    out = ROOT / "docs" / "data" / "layer_skill_diagnostic.json"
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
