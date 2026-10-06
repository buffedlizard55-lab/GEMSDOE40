#!/usr/bin/env python3
"""H41 — multi-scale-stable shallow-contact Euler depth-cluster emission.

Pre-registration (frozen before any scoring; see ``docs/research/h41-preregistration-20261006.md``):

1. **Euler deconvolution** (Reid, Allsop, Granser, Millett & Somerton 1990, *Geophysics*
   55(1) 80-91, doi:10.1190/1.1442774) at the fault-like **contact** index SI = 0 on the
   organizer's ``rtp`` (band 2), ``tmi`` (band 14) and ``iso_grav_anom`` (band 13) layers,
   moving-window least squares with the repository's existing gating (analytic-signal
   percentile, conditioning, source-in-window, relative depth standard error, depth band).
   Outputs a cloud of depth-labelled solution points, never an edge map.
2. **Depth-cluster KDE** (:mod:`gemsdoe40.kde`): weighted bilinear splat + Gaussian KDE,
   weights = shallowness x solution quality x local tightness x mutual depth consistency,
   plus a magnetic/gravity concordance term where both families agree in position *and*
   depth.  Additional H41 restriction: only solutions shallower than ``--max-depth`` enter
   the KDE (the brief's "tight clusters of shallow, mutually-consistent solutions").
3. **Multi-scale stability (new in H41).**  Each window length produces its own KDE.
   The emitted score is the geometric mean across window lengths, so a cell must be
   supported by *every* window scale rather than by one lucky window; the per-window
   count is recorded.  This is what removes single-window artefacts that a KDE of a
   merged cloud cannot distinguish.
4. **Metric-optimal emission** (:mod:`gemsdoe40.emission`): value-ranked Poisson-disk
   thinning at 3 px spacing (>= the metric's own support R = 300 m, so no two dots can
   claim the same truth pixel) and a mass budget chosen from the recalibrated
   live-score instrument, subject to the exact marginal rule of the published metric.
5. **Audit** before publishing: exact format, credit-bar check, blocked off-catalogue
   proxy comparison against the owner's two best scored dot sets, and a
   hash/Pearson/Jaccard near-duplicate test against every staged prior raster.

Writes ``docs/downloads/``-ready GeoTIFFs into ``work/`` first; nothing is published and
nothing is uploaded.  ``--publish`` copies the audited bytes into ``docs/downloads/``.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import pickle
import sys
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import rasterio
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe40 import GRID_HEIGHT, GRID_WIDTH  # noqa: E402
from gemsdoe40.emission import (  # noqa: E402
    calibrated_live_score,
    value_ranked_thinning,
)
from gemsdoe40.euler import EulerCloud, deconvolve, merge_clouds  # noqa: E402
from gemsdoe40.grid import footprint_from_sample, read_band, read_labels, sha256  # noqa: E402
from gemsdoe40.grid import write_submission  # noqa: E402
from gemsdoe40.discriminant import auc, blocked_folds, fit_logistic, predict_logistic  # noqa: E402
from gemsdoe40.kde import concordance_boost, solution_weights, splat_kde  # noqa: E402
from gemsdoe40.measure import credit_components, kernel_field  # noqa: E402

DATA = ROOT / "data"
WORK = ROOT / "work"
DOWNLOADS = ROOT / "docs" / "downloads"
REF_PRIOR = ROOT / "ref" / "prior"

FIELDS = (("rtp", 2, "magnetic"), ("tmi", 14, "magnetic"), ("iso_grav_anom", 13, "gravity"))
WINDOWS = (10, 16, 24)
STRIDE = 4
DEPTH_MAX_M = 1500.0
SPACING_PX = 3.0
CATALOGUE_FLANK_PX = 2.0
MASS_GRID = (10_000, 20_000, 30_000, 40_000, 60_000, 90_000, 120_000)
instrument = None  # filled from docs/data/instrument-recalibration.json


def log(msg: str) -> None:
    print(msg, flush=True)


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def instrument_score(mass: int, w: float) -> float:
    if instrument is None:
        return calibrated_live_score(mass, w)
    K, a = instrument["K"], instrument["a"]
    tp = K * (1.0 - np.exp(-a * mass * w / K))
    return float(tp / (0.2 * tp + 0.2 * mass + 0.8 * K))


# --------------------------------------------------------------------------- #
# 1 + 2 + 3: Euler cloud -> per-window KDE -> scale-stable field
# --------------------------------------------------------------------------- #

def euler_clouds(cache: Path | None = None) -> dict[int, list[EulerCloud]]:
    """Deconvolve every field at every window; return clouds grouped by window."""
    key = {"fields": FIELDS, "windows": list(WINDOWS), "stride": STRIDE,
           "si": 0.0, "features_sha": sha256(DATA / "training_features.tif")}
    if cache is not None and cache.exists():
        blob = pickle.loads(cache.read_bytes())
        if blob.get("key") == key:
            log(f"  Euler clouds restored from {cache}")
            return blob["by_window"]
    footprint = footprint_from_sample(DATA / "sample_submission.tif")
    by_window: dict[int, list[EulerCloud]] = {w: [] for w in WINDOWS}
    for name, band, family in FIELDS:
        arr, valid = read_band(DATA / "training_features.tif", band)
        valid = valid & footprint
        for win in WINDOWS:
            t0 = time.time()
            cloud = deconvolve(
                arr, valid, field_name=name, structural_index=0.0,
                window_px=win, stride_px=STRIDE, analytic_percentile=72.0,
                max_rel_se=0.22, min_depth_m=80.0, max_depth_m=2200.0,
            )
            cloud.stats["family"] = family
            by_window[win].append(cloud)
            log(f"  Euler {name:14s} win={win:2d} SI=0 -> {len(cloud):7,d} solutions "
                f"({time.time() - t0:5.1f} s, gated {cloud.stats.get('gated_windows', 0):,})")
        del arr, valid
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_bytes(pickle.dumps({"key": key, "by_window": by_window}))
        log(f"  cached Euler clouds -> {cache}")
    return by_window


def shallow(cloud: EulerCloud, depth_max: float) -> EulerCloud:
    keep = cloud.depth_m <= depth_max
    return EulerCloud(
        cloud.field_name, cloud.structural_index, cloud.window_px, cloud.stride_px,
        cloud.row[keep], cloud.col[keep], cloud.depth_m[keep], cloud.rel_se[keep],
        cloud.cond[keep], cloud.analytic[keep],
        stats={**cloud.stats, "kept_shallow": int(keep.sum()), "depth_max_m": depth_max},
    )


def depth_map(clouds: list[EulerCloud], shape, *, sigma: float, depth_max: float) -> np.ndarray:
    """Weighted mean solution depth per cell (0 where no solution supports the cell)."""
    num = np.zeros(shape, dtype=np.float64)
    den = np.zeros(shape, dtype=np.float64)
    for c in clouds:
        cs = shallow(c, depth_max)
        w = solution_weights(cs)
        num += splat_kde(cs, w * np.clip(cs.depth_m, 0, None), shape, sigma_px=sigma)
        den += splat_kde(cs, w, shape, sigma_px=sigma)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(den > 0, num / np.maximum(den, 1e-12), 0.0)
    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)


def window_field(clouds: list[EulerCloud], shape, *, sigma: float, concord: float,
                 depth_max: float) -> tuple[np.ndarray, dict]:
    """One window's depth-cluster KDE: magnetic + gravity splats + concordance."""
    total = np.zeros(shape, dtype=np.float64)
    per = {}
    for c in clouds:
        cs = shallow(c, depth_max)
        w = solution_weights(cs)
        f = splat_kde(cs, w, shape, sigma_px=sigma)
        total += f
        per[f"{c.field_name}:{c.window_px}"] = {
            "solutions": int(len(c)), "kept_shallow": int(len(cs)),
            "weight_sum": float(w.sum()), "field_max": float(f.max()),
        }
    mag = merge_clouds([shallow(c, depth_max) for c in clouds
                        if c.field_name in ("rtp", "tmi")])
    grav = merge_clouds([shallow(c, depth_max) for c in clouds
                         if c.field_name == "iso_grav_anom"])
    conc = concordance_boost(mag, grav, shape, sigma_px=max(1.2, sigma - 0.4))
    total = total + concord * conc
    return total, {"per_cloud": per, "concordance_max": float(conc.max())}


def scale_stable(fields: dict[int, np.ndarray]) -> tuple[np.ndarray, dict, dict]:
    """Geometric mean of the per-window normalised KDEs (0 where any window is 0)."""
    norm = {}
    for win, f in fields.items():
        peak = float(f.max())
        norm[win] = (f / peak) if peak > 0 else f
    stack = np.stack([norm[w] for w in sorted(norm)])
    supported = np.ones(stack.shape[1:], dtype=bool)
    for s in stack:
        supported &= s > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        geo = np.exp(np.mean(np.log(np.where(stack > 0, stack, np.nan)), axis=0))
    geo = np.where(supported, geo, 0.0)
    counts = (stack > 0).sum(axis=0)
    return geo.astype(np.float32), {
        "windows": sorted(norm),
        "cells_supported_all_windows": int(supported.sum()),
        "cells_supported_any": int((counts > 0).sum()),
        "per_window_max": {str(w): float(norm[w].max()) for w in norm},
        "geometric_mean_max": float(geo.max()),
    }, norm


# --------------------------------------------------------------------------- #
# 3b: spatially blocked discriminant trained on the off-catalogue public fault map
# --------------------------------------------------------------------------- #

FEATURE_BANDS = (("rtp", 2), ("tmi_hg", 3), ("geod_2ndinv", 4), ("iso_grav_anom_slope", 5),
                 ("tc", 6), ("iso_grav_anom_vg", 11), ("det_elev", 12), ("iso_grav_anom", 13),
                 ("tmi", 14), ("depth_to_base_surf", 15))
#: Surface-expression features derived from the detrended-elevation band.  The 300 m
#: metric credit and the observed leaderboard behaviour of structure-aligned dot sets
#: both point at lineaments, so the discriminant gets curvature at two smoothing scales
#: in addition to the raw band.
DEM_SCALES_PX = (1.0, 3.0)
RATIO_BANDS = (("vg_ratio", 9, 11), ("hg_ratio", 3, 5))
FEATURE_NAMES = (["euler_w%d" % w for w in WINDOWS] + ["euler_scale_stable"]
                 + ["euler_depth_w%d" % w for w in WINDOWS]
                 + ["grad_rtp", "grad_grav"]
                 + ["dem_grad"] + [f"dem_lap_s{int(s)}" for s in DEM_SCALES_PX]
                 + [n for n, _ in FEATURE_BANDS]
                 + [n for n, _, _ in RATIO_BANDS] + ["rtp_lap_s2", "tc_z15", "dist_catalogue"])
DILATE_POS_PX = 1.0
CLEAN_NEG_PX = 3.0
MAX_TRAIN_POS = 60_000
MAX_TRAIN_NEG = 240_000


def build_feature_matrix(support_idx: tuple[np.ndarray, np.ndarray],
                         windows_norm: dict[int, np.ndarray],
                         scale_field: np.ndarray,
                         truth: np.ndarray,
                         depth_maps: dict[int, np.ndarray] | None = None
                         ) -> tuple[np.ndarray, np.ndarray, list[str], int]:
    """Feature rows at the raster cells listed by ``support_idx`` (row, col arrays).

    Every geophysical band is read once and reduced to those cells, so the peak
    footprint is one band plus its local gradients rather than the whole stack.
    """
    rows, cols = support_idx
    cols_list: list[np.ndarray] = []
    for win in WINDOWS:
        cols_list.append(np.log1p(np.asarray(windows_norm[win], dtype=np.float32)[rows, cols] * 10.0))
    cols_list.append(np.log1p(np.asarray(scale_field, dtype=np.float32)[rows, cols] * 10.0))
    for win in WINDOWS:
        dm = depth_maps[win] if depth_maps else np.zeros_like(scale_field)
        cols_list.append(np.asarray(dm, dtype=np.float32)[rows, cols].astype(np.float64) / 1000.0)
    for band_name, band in (("rtp", 2), ("iso_grav_anom", 13)):
        arr, valid = read_band(DATA / "training_features.tif", band)
        # zero-fill invalid cells *before* differencing so no NaN propagates;
        # the sample cells are all valid (they lie inside the footprint).
        arr = np.where(valid, arr, 0.0).astype(np.float32)
        gr, gc = np.gradient(arr)
        g = np.hypot(gr, gc).astype(np.float32)
        cols_list.append(np.log1p(g[rows, cols]))
        del arr, g, gr, gc, valid
    # detrended-elevation surface expression: gradient magnitude and two-scale curvature
    dem, dem_valid = read_band(DATA / "training_features.tif", 12)
    dem = np.where(dem_valid, dem, 0.0).astype(np.float32)
    gr, gc = np.gradient(dem)
    cols_list.append(np.log1p(np.hypot(gr, gc)[rows, cols]).astype(np.float64))
    del gr, gc
    for scale in DEM_SCALES_PX:
        sm = ndimage.gaussian_filter(dem, scale) if scale > 0 else dem
        lap = ndimage.laplace(sm)
        cols_list.append(np.log1p(np.abs(lap)[rows, cols]).astype(np.float64))
        del sm, lap
    del dem, dem_valid
    for name, band in FEATURE_BANDS:
        arr, valid = read_band(DATA / "training_features.tif", band)
        cols_list.append(np.where(valid, arr, 0.0)[rows, cols].astype(np.float64))
        # missing cells inside a support cell are imputed with 0 after z-scoring
        del arr, valid
    # H41b additions: cross-field vertical-gradient ratios (the H43 hypothesis), the
    # short-scale curvature of the reduced-to-pole field, a local radiometric z-score and
    # the distance to the provided catalogue (fault networks are spatially clustered:
    # surrogate fault pixels sit at mean kernel 0.149 from the catalogue versus 0.050 for
    # random cells, so this is real, non-circular signal about *where* faults are).
    for _, num_band, den_band in RATIO_BANDS:
        num, nvalid = read_band(DATA / "training_features.tif", num_band)
        den, dvalid = read_band(DATA / "training_features.tif", den_band)
        num = np.where(nvalid, num, 0.0).astype(np.float32)
        den = np.where(dvalid, den, 0.0).astype(np.float32)
        ratio = np.abs(num) / (np.abs(den) + 1e-3)
        cols_list.append(np.log1p(ratio[rows, cols]).astype(np.float64))
        del num, den, nvalid, dvalid, ratio
    rtp, rvalid = read_band(DATA / "training_features.tif", 2)
    rtp = np.where(rvalid, rtp, 0.0).astype(np.float32)
    cols_list.append(np.log1p(np.abs(ndimage.laplace(ndimage.gaussian_filter(rtp, 2.0)))[rows, cols])
                     .astype(np.float64))
    del rtp, rvalid
    tc, tvalid = read_band(DATA / "training_features.tif", 6)
    tc = np.where(tvalid, tc, 0.0).astype(np.float32)
    mean = ndimage.uniform_filter(tc, 15)
    sq = ndimage.uniform_filter(tc * tc, 15)
    std = np.sqrt(np.maximum(sq - mean * mean, 0.0)) + 1e-3
    cols_list.append(((tc - mean) / std)[rows, cols].astype(np.float64))
    del tc, tvalid, mean, sq, std
    dcat = ndimage.distance_transform_edt(~read_labels(DATA / "labels.tif"))
    cols_list.append(np.log1p(np.minimum(dcat, 50.0))[rows, cols].astype(np.float64))
    del dcat
    x = np.column_stack(cols_list).astype(np.float64)
    n_bad = int((~np.isfinite(x)).sum())
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    y = truth[rows, cols].astype(np.float64)
    return x, y, list(FEATURE_NAMES), n_bad


def oof_discriminant(support: np.ndarray, truth: np.ndarray, windows_norm, scale_field,
                     *, depth_maps: dict[int, np.ndarray] | None = None,
                     rows_blocks: int = 5, cols_blocks: int = 6, seed: int = 20261006,
                     l2: float = 1.0) -> tuple[np.ndarray, dict]:
    """Out-of-fold fault probability on the Euler support, with a blocked audit."""
    from scipy import ndimage

    rng = np.random.default_rng(seed)
    rr, cc = np.nonzero(support)
    n_support = rr.size
    x_all, y_all, names, n_bad = build_feature_matrix((rr, cc), windows_norm, scale_field,
                                                     truth, depth_maps)

    fold_id, n_rows, n_cols = blocked_folds(support.shape, rows_blocks, cols_blocks)
    f_all = fold_id[rr, cc]

    pos_band = ndimage.binary_dilation(truth, iterations=int(np.ceil(DILATE_POS_PX)))
    dist = ndimage.distance_transform_edt(~truth)
    credit = np.maximum(1.0 - dist / 3.0, 0.0)
    y_credit = credit[rr, cc].astype(np.float64)
    y_pos = pos_band[rr, cc] & (dist[rr, cc] <= DILATE_POS_PX)
    y_clean = (~pos_band[rr, cc]) & (dist[rr, cc] > CLEAN_NEG_PX)

    pos_idx = np.nonzero(y_pos)[0]
    neg_idx = np.nonzero(y_clean)[0]
    if pos_idx.size > MAX_TRAIN_POS:
        pos_idx = rng.choice(pos_idx, MAX_TRAIN_POS, replace=False)
    if neg_idx.size > MAX_TRAIN_NEG:
        neg_idx = rng.choice(neg_idx, MAX_TRAIN_NEG, replace=False)

    prob = np.zeros(n_support, dtype=np.float32)
    folds: list[dict] = []
    for fold in range(n_rows * n_cols):
        train = np.zeros(n_support, dtype=bool)
        train[pos_idx] = True
        train[neg_idx] = True
        in_fold = f_all == fold
        train &= ~in_fold
        if train.sum() < 50 or not in_fold.any():
            prob[in_fold] = float(y_pos[np.nonzero(in_fold)[0]].mean()) if in_fold.any() else 0.0
            continue
        sw = np.where(y_pos[train], 1.0, 1.0)
        soft = np.where(y_pos[train], np.maximum(y_credit[train], 0.34), 1.0)
        model = fit_logistic(x_all[train], y_pos[train].astype(float), l2=l2,
                             sample_weight=soft * sw)
        held = np.nonzero(in_fold)[0]
        prob[held] = predict_logistic(model, x_all[held]).astype(np.float32)
        # honest leave-this-block-out AUC on the held block's labelled sample
        lab = np.isin(held, np.concatenate([pos_idx, neg_idx]))
        if lab.sum() > 20:
            folds.append({"fold": int(fold), "n_test": int(lab.sum()),
                          "auc": auc(y_pos[held[lab]], prob[held[lab]]),
                          "pos_rate": float(y_pos[held[lab]].mean())})
    lab_all = np.zeros(n_support, dtype=bool)
    lab_all[pos_idx] = True
    lab_all[neg_idx] = True
    overall = auc(y_pos[lab_all], prob[lab_all])
    field = np.zeros(support.shape, dtype=np.float32)
    field[rr, cc] = prob
    report = {
        "n_support_px": int(n_support),
        "n_train_pos_available": int(y_pos.sum()),
        "n_train_neg_available": int(y_clean.sum()),
        "n_train_pos_sampled": int(pos_idx.size),
        "n_train_neg_sampled": int(neg_idx.size),
        "features": names,
        "nonfinite_feature_cells_replaced": n_bad,
        "oof_auc_all_blocks": overall,
        "oof_auc_mean_of_blocks": float(np.mean([f["auc"] for f in folds if f["auc"]])) if folds else None,
        "blocks_with_auc": [f for f in folds if f["auc"]],
        "positive_label": (f"proxy fault pixels dilated by {DILATE_POS_PX} px with per-sample "
                           f"weight = triangular kernel credit (1.0/0.67/0.34); negatives are "
                           f">{CLEAN_NEG_PX} px from any proxy fault; the 1-3 px band is dropped"),
        "proxy_warning": ("labels are SGMC faults absent from the provided catalogue: an "
                          "independent public stand-in for the unmapped-fault class, NOT the "
                          "organizer's hidden truth"),
    }
    return field, report


# --------------------------------------------------------------------------- #
# 4: emission
# --------------------------------------------------------------------------- #

def emission_curve(score: np.ndarray, support: np.ndarray, truth: np.ndarray,
                   footprint: np.ndarray) -> list[dict]:
    """Streaming value-ranked thinning: metric components at each mass snapshot."""
    from gemsdoe40.emission import thinning_curve_snapshots

    rows: list[dict] = []

    def snapshot(mass: int, emit: np.ndarray) -> None:
        comp = credit_components(emit, truth, valid=footprint)
        n = int(mass)
        w = comp["pred_credit"] / n if n else 0.0
        rows.append({
            "mass": n, "spacing_px": SPACING_PX,
            "proxy_tp": comp["tp"], "proxy_dti": comp["dti"],
            "credit_per_dot_w": w, "credit_per_mass_u": comp["credit_per_mass"],
            "instrument_pred": instrument_score(n, w),
        })
        log(f"    mass {n:8,d}  proxyDTI {comp['dti']:.4f}  w {w:.4f}  "
            f"instrument {rows[-1]['instrument_pred']:.4f}")

    thinning_curve_snapshots(score, support, spacing_px=SPACING_PX,
                             breakpoints=list(MASS_GRID), on_snapshot=snapshot)
    return rows


def blocked_compare(emit: np.ndarray, footprint: np.ndarray, truth: np.ndarray,
                    catalogue: np.ndarray, comparators: dict[str, np.ndarray],
                    *, rows_blocks: int = 4, cols_blocks: int = 6, guard: int = 3) -> dict:
    """Spatially blocked DTI of the candidate against each comparator.

    Blocks are ``rows_blocks x cols_blocks``; a ``guard``-pixel border is removed from
    every block so a prediction outside a block cannot earn credit inside it.  Every
    raster is scored on the identical instrument.
    """
    h, w = footprint.shape
    out: dict[str, dict] = {}
    bh, bw = h // rows_blocks, w // cols_blocks
    cand = {"H41-candidate": np.asarray(emit, dtype=np.float32),
            **{k: np.asarray(v, dtype=np.float32) for k, v in comparators.items()}}
    per_block = {name: [] for name in cand}
    kept = []
    for i in range(rows_blocks):
        for j in range(cols_blocks):
            r0, r1 = i * bh + guard, (i + 1) * bh - guard
            c0, c1 = j * bw + guard, (j + 1) * bw - guard
            if r1 <= r0 or c1 <= c0:
                continue
            valid = np.zeros_like(footprint)
            valid[r0:r1, c0:c1] = True
            valid &= footprint
            t = truth & valid & ~catalogue
            if int(t.sum()) < 25:
                continue
            kept.append((i, j, int(t.sum())))
            for name, pred in cand.items():
                comp = credit_components(pred * valid, t, valid=valid)
                gc.collect()
                per_block[name].append({"block": f"r{i}c{j}", "dti": comp["dti"],
                                        "tp": comp["tp"], "truth_px": int(t.sum())})
    for name, bl in per_block.items():
        vals = np.array([b["dti"] for b in bl], float)
        out[name] = {"per_block": bl, "mean_dti": float(vals.mean()) if vals.size else 0.0,
                     "median_dti": float(np.median(vals)) if vals.size else 0.0}
    wins = sum(1 for b, c in zip(per_block["H41-candidate"],
                                 zip(*(per_block[n] for n in per_block if n != "H41-candidate")))
               if b["dti"] > max(x["dti"] for x in c))
    out["summary"] = {
        "blocks_scored": len(kept),
        "blocks": [f"r{i}c{j}:{n}" for i, j, n in kept],
        "strict_wins_of_candidate": int(wins),
        "comparators": [n for n in per_block if n != "H41-candidate"],
    }
    return out


# --------------------------------------------------------------------------- #
# 5: novelty
# --------------------------------------------------------------------------- #

def novelty(emit: np.ndarray, footprint: np.ndarray, paths: list[Path], top_k: int) -> dict:
    rows = []
    for p in sorted(paths):
        try:
            with rasterio.open(p) as ds:
                a = ds.read(1).astype(np.float32)
        except Exception as exc:
            rows.append({"file": p.name, "error": str(exc)})
            continue
        if a.shape != emit.shape:
            rows.append({"file": p.name, "error": f"shape {a.shape}"})
            continue
        prior = np.where(np.isfinite(a), a, 0.0)
        prior_bin = (prior > 0) & footprint
        fv, pv = emit[footprint].astype(np.float32), prior[footprint]
        pearson = float(np.corrcoef(fv, pv)[0, 1]) if fv.std() > 0 and pv.std() > 0 else 0.0
        inter = int((emit.astype(bool) & prior_bin).sum())
        union = int((emit.astype(bool) | prior_bin).sum())
        jac = inter / union if union else 0.0
        top_prior = np.zeros_like(emit, bool)
        if top_k:
            idx = np.argsort(-prior.ravel(), kind="stable")[:top_k]
            top_prior.ravel()[idx] = True
            top_prior &= footprint
        top_inter = int((emit.astype(bool) & top_prior).sum())
        top_union = int((emit.astype(bool) | top_prior).sum())
        rows.append({"file": p.name, "sha256": sha256(p),
                     "prior_positive_px": int(prior_bin.sum()),
                     "pearson": pearson, "support_jaccard": jac,
                     "top_budget_jaccard": (top_inter / top_union) if top_union else 0.0})
    scored = [r for r in rows if "pearson" in r]
    worst_p = max(scored, key=lambda r: abs(r["pearson"])) if scored else None
    worst_j = max(scored, key=lambda r: r["support_jaccard"]) if scored else None
    worst_t = max(scored, key=lambda r: r["top_budget_jaccard"]) if scored else None
    return {
        "n_priors": len(rows),
        "max_abs_pearson": worst_p["pearson"] if worst_p else None,
        "max_abs_pearson_file": worst_p["file"] if worst_p else None,
        "max_support_jaccard": worst_j["support_jaccard"] if worst_j else None,
        "max_support_jaccard_file": worst_j["file"] if worst_j else None,
        "max_top_budget_jaccard": worst_t["top_budget_jaccard"] if worst_t else None,
        "max_top_budget_jaccard_file": worst_t["file"] if worst_t else None,
        "is_new": bool(worst_p and abs(worst_p["pearson"]) < 0.85
                       and worst_j["support_jaccard"] < 0.35
                       and worst_t["top_budget_jaccard"] < 0.35),
        "rows": rows,
    }


def main() -> int:
    global instrument
    ap = argparse.ArgumentParser()
    ap.add_argument("--publish", action="store_true")
    ap.add_argument("--rank", choices=("euler", "gated"), default="gated")
    ap.add_argument("--flank", type=float, default=2.0,
                    help="catalogue exclusion flank in px (2 px = 200 m, the measured "
                         "minimum distance of the best live-scoring prior artifact)")
    ap.add_argument("--slug", default=None)
    ap.add_argument("--cache", type=Path, default=ROOT / "work" / "h41_euler_clouds.pkl")
    args = ap.parse_args()

    global CATALOGUE_FLANK_PX
    CATALOGUE_FLANK_PX = float(args.flank)

    t0 = time.time()
    WORK.mkdir(exist_ok=True)
    receipt: dict = {"experiment": "H41", "started_utc": now(),
                     "method": "multi-scale-stable shallow-contact Euler depth-cluster emission",
                     "windows_px": list(WINDOWS), "stride_px": STRIDE,
                     "depth_max_m": DEPTH_MAX_M, "spacing_px": SPACING_PX,
                     "catalogue_flank_px": CATALOGUE_FLANK_PX}

    recal = ROOT / "docs" / "data" / "instrument-recalibration.json"
    if recal.exists():
        blob = json.loads(recal.read_text(encoding="utf-8"))
        instrument = {"K": blob["refit"]["K"], "a": blob["refit"]["a"],
                      "loo_spearman": blob["refit"]["loo_spearman"],
                      "source": "docs/data/instrument-recalibration.json (refit on measured anchors)"}
        log(f"instrument: K={instrument['K']:.1f} a={instrument['a']:.4f} "
            f"LOO Spearman={instrument['loo_spearman']}")
    receipt["instrument"] = instrument

    footprint = footprint_from_sample(DATA / "sample_submission.tif")
    catalogue = read_labels(DATA / "labels.tif")
    with rasterio.open(DATA / "external" / "derived_sgmc_faults_100m_u8.tif") as ds:
        sgmc = ds.read(1) > 0
    truth_proxy = sgmc & ~catalogue & footprint
    log(f"grid {footprint.shape} footprint {int(footprint.sum()):,} "
        f"catalogue {int(catalogue.sum()):,} proxy truth {int(truth_proxy.sum()):,}")
    receipt["grid"] = {"shape": list(footprint.shape), "footprint_px": int(footprint.sum()),
                       "catalogue_px": int(catalogue.sum()),
                       "proxy_truth_px": int(truth_proxy.sum())}

    log("1. Euler deconvolution (Reid et al. 1990), SI = 0 fault-like contact")
    by_window = euler_clouds(args.cache)
    receipt["euler"] = {
        "per_window": {str(w): {"solutions": int(sum(len(c) for c in cs)),
                                "per_field": {c.field_name: int(len(c)) for c in cs}}
                       for w, cs in by_window.items()},
        "total_solutions": int(sum(len(c) for cs in by_window.values() for c in cs)),
    }
    for w, cs in by_window.items():
        log(f"   window {w:2d} px: {sum(len(c) for c in cs):,} solutions")

    log("2/3. depth-cluster KDE per window + multi-scale geometric mean")
    fields, stats, depths = {}, {}, {}
    for w, cs in by_window.items():
        f, s = window_field(cs, footprint.shape, sigma=1.7, concord=1.4, depth_max=DEPTH_MAX_M)
        fields[w] = np.where(footprint, f, 0.0)
        depths[w] = depth_map(cs, footprint.shape, sigma=1.7, depth_max=DEPTH_MAX_M)
        stats[str(w)] = s
        log(f"   window {w:2d}: KDE max {fields[w].max():.4f} "
            f"positive {int((fields[w] > 0).sum()):,}")
    field, stability, norm_fields = scale_stable(fields)
    receipt["kde"] = {"per_window": stats, "stability": stability}
    log(f"   scale-stable field: max {field.max():.4f} "
        f"supported-all-windows {stability['cells_supported_all_windows']:,}")

    support = footprint & ~catalogue
    if CATALOGUE_FLANK_PX > 0:
        dcat = ndimage.distance_transform_edt(~catalogue)
        support &= dcat > CATALOGUE_FLANK_PX
    # The manifest requirement of this arm: every emitted dot lies on a cell with
    # Euler depth-cluster support at every window scale.
    support &= field > 0
    receipt["support_px"] = int(support.sum())
    log(f"   Euler-gated support {int(support.sum()):,} px")

    rank_field = field
    if args.rank == "gated":
        log("3b. spatially blocked discriminant on off-catalogue public faults")
        rank_field, disc = oof_discriminant(support, truth_proxy, norm_fields, field,
                                            depth_maps=depths)
        receipt["discriminant"] = disc
        log(f"   OOF AUC {disc['oof_auc_all_blocks']:.4f} "
            f"(mean of blocks {disc['oof_auc_mean_of_blocks']:.4f}, "
            f"{len(disc['blocks_with_auc'])} blocks)")
        rank_field = np.where(support, rank_field, 0.0).astype(np.float32)

    log("4. metric-optimal emission")
    curve = emission_curve(rank_field, support, truth_proxy, footprint)
    receipt["curve"] = curve
    if not curve:
        raise SystemExit("no emission candidates: ranking field is empty or non-finite")
    best = max(curve, key=lambda r: (r["instrument_pred"], r["proxy_dti"]))
    emit = value_ranked_thinning(rank_field, support, spacing_px=SPACING_PX,
                                 max_mass=int(best["mass"]))
    mass = int(emit.sum())
    log(f"   chosen: mass {mass:,} at {SPACING_PX} px spacing "
        f"(instrument {best['instrument_pred']:.4f}, proxy DTI {best['proxy_dti']:.4f})")
    receipt["chosen"] = {"mass": mass, "spacing_px": SPACING_PX, **best}

    # credit-bar check on the engine's own emitted dots (exact, off-catalogue domain).
    kt = np.maximum(1.0 - ndimage.distance_transform_edt(~truth_proxy) / 3.0, 0.0)
    w_emit = float(kt[emit > 0].mean()) if mass else 0.0
    receipt["emitted"] = {"mass": mass, "credit_per_dot_w": w_emit,
                          "instrument_pred": instrument_score(mass, w_emit),
                          "rank": args.rank}
    if args.rank == "gated":
        pref = float(((emit > 0) & (kt > 0)).sum()) / max(mass, 1)
        receipt["emitted"]["fraction_of_dots_within_3px_of_proxy"] = pref

    log("5. audit")
    comp = credit_components(emit, truth_proxy, valid=footprint)
    cat_comp = credit_components(emit, catalogue, valid=footprint)
    receipt["audit"] = {
        "proxy_dti": comp["dti"], "proxy_tp": comp["tp"], "proxy_fp": comp["fp"],
        "proxy_fn": comp["fn"], "proxy_truth_px": comp["n_truth"],
        "credit_per_mass_u": comp["credit_per_mass"],
        "catalogue_dti": cat_comp["dti"], "dots_on_catalogue": int(cat_comp["tp"]),
    }
    log(f"   proxy DTI {comp['dti']:.4f} (T {comp['tp']:.0f}, F {comp['fp']:.0f}, "
        f"FN {comp['fn']:.0f}); catalogue DTI {cat_comp['dti']:.4f}, "
        f"dots-on-catalogue {int(cat_comp['tp'])}")

    comparators: dict[str, np.ndarray] = {}
    for label, name in (("H33-B2", "gemsdoe32-h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros.tif"),
                        ("H27-4", "gems28-h27-4-r1-solo-d2-8-20261003-8acb75e1f2cc-allfinite.tif"),
                        ("H40-E", "gemsdoe39-h40-e-disc-h40e-30k-zeros.tif")):
        p = REF_PRIOR / name
        if p.exists():
            with rasterio.open(p) as ds:
                a = ds.read(1)
            comparators[label] = ((np.isfinite(a) & (a > 0)) & footprint).astype(np.float32)
            cc = credit_components(comparators[label], truth_proxy, valid=footprint)
            comparators[label + "_proxy_dti"] = np.array(cc["dti"])  # keep the number
    receipt["comparators"] = {k: (float(v) if v.ndim == 0 else int(v.sum()))
                              for k, v in comparators.items() if k.endswith("_proxy_dti")}
    blocked = blocked_compare(emit, footprint, truth_proxy, catalogue,
                              {k: v for k, v in comparators.items() if not k.endswith("_proxy_dti")})
    receipt["blocked"] = blocked
    log(f"   blocked: candidate wins {blocked['summary']['strict_wins_of_candidate']} of "
        f"{blocked['summary']['blocks_scored']} truth-bearing blocks vs "
        f"{blocked['summary']['comparators']}")

    nov = novelty(emit, footprint, sorted(REF_PRIOR.glob("*.tif")), top_k=mass)
    receipt["novelty"] = {k: v for k, v in nov.items() if k != "rows"}
    log(f"   novelty: max|r| {nov['max_abs_pearson']:.4f}, max Jaccard "
        f"{nov['max_support_jaccard']:.4f}, top-budget Jaccard "
        f"{nov['max_top_budget_jaccard']:.4f} -> is_new={nov['is_new']}")

    # write both portal conventions into work/
    slug = args.slug or f"gemsdoe40-h41-msst-euler-si0-depthcluster-{now()}"
    digest = hashlib.sha256(np.ascontiguousarray(emit.astype(np.float32)).tobytes()).hexdigest()
    name = f"{slug}-{digest[:8]}"
    zeros = WORK / f"{name}-zeros.tif"
    nan = WORK / f"{name}-nan.tif"
    write_submission(emit, DATA / "sample_submission.tif", zeros, footprint, outside="zero")
    write_submission(emit, DATA / "sample_submission.tif", nan, footprint, outside="nan")
    with zipfile.ZipFile(WORK / f"{name}-zeros.zip", "w", zipfile.ZIP_DEFLATED) as z:
        z.write(zeros, zeros.name)
    receipt["files"] = {
        "slug": name,
        "zeros": {"path": str(zeros.relative_to(ROOT)), "sha256": sha256(zeros),
                  "bytes": zeros.stat().st_size},
        "nan": {"path": str(nan.relative_to(ROOT)), "sha256": sha256(nan),
                "bytes": nan.stat().st_size},
    }
    receipt["finished_utc"] = now()
    receipt["seconds"] = round(time.time() - t0, 1)

    out = ROOT / "docs" / "data" / "h41-generation.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    (WORK / "h41-uniqueness-full.json").write_text(json.dumps(nov, indent=2) + "\n",
                                                   encoding="utf-8")
    log(f"[out] {out}")
    log(json.dumps({k: v for k, v in receipt.items()
                    if k not in ("curve", "kde", "euler", "novelty", "blocked")}, indent=2))
    if args.publish:
        DOWNLOADS.mkdir(parents=True, exist_ok=True)
        for p in (zeros, nan, WORK / f"{name}-zeros.zip"):
            target = DOWNLOADS / p.name
            target.write_bytes(p.read_bytes())
            log(f"published {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
