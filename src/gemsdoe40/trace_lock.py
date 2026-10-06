"""H8: trace-locked Euler depth-consensus emission.

The H4 experiment proved the solver physics produces a real, depth-labelled solution cloud but a
diffuse KDE emission (865,145 positive cells) that is not spatially selective
(`docs/reports/h4a-negative-result-20261006.md`, AUC 0.534 vs catalogue).  Every artifact in the
corpus that scored >= 0.24 live is trace/dot-like rather than blob-like.  H8 keeps the frozen H4
cloud and changes only the *emission geometry*:

1. per-family weighted KDE of the cluster-weighted cloud (frozen H4 settings);
2. combined consensus field ``C = (M + G + sqrt(M G))/3`` (frozen H4 rule);
3. orientation from the structure tensor of ``C``;
4. orientation-NMS ridge crests (1-px thinning, plateau-safe one-sided tie rule);
5. along-strike depth consensus: a crest passes only if its oriented corridor
   (+/-8 px along strike, +/-3 px across) contains >= 3 weighted solutions with corridor depth
   standard deviation <= 445 m (normal-equivalent of MAD <= 300 m);
6. continuous output ``gaussian(C * gamma, sigma=1)`` normalized to [0,1], exact known pixels
   zeroed, NaN outside the sample footprint.

All constants are frozen in ``docs/research/h8-preregistration-20261006.md``.  No proxy truth,
holdout label, prior prediction or leaderboard value enters construction.
"""
from __future__ import annotations

import csv
import gzip
import io
from pathlib import Path

import numpy as np
from scipy import ndimage

from .contact_euler import bilinear_splat

CLOUD_SHA256 = "6bed30b224981ef1064256d696c585b76f438011bbae6c0a16c00e6e6f5922ba"
N_RECORDS = 46_656
FAMILIES = {"tmi": 0, "iso_grav_anom": 1}

# Frozen emission constants (preregistration sections 4-7).
KDE_SIGMA_PX = 2.0
KDE_TRUNCATE = 4.0
FAMILY_PERCENTILE = 99.5
TENSOR_SIGMA_INNER_PX = 1.0
TENSOR_SIGMA_INTEGRATION_PX = 3.0
ISOTROPY_MIN = 0.15
CORRIDOR_ALONG_HALF_PX = 8
CORRIDOR_ACROSS_HALF_PX = 3
MIN_CORRIDOR_SOLUTIONS = 3
CORRIDOR_SD_MAX_M = 445.0
OUTPUT_SIGMA_PX = 1.0
OUTPUT_TRUNCATE = 3.0

#: gradient direction per quantized bin, as (d_row, d_col); orientation is mod 180 degrees.
BIN_DIRECTIONS = {0: (0, 1), 1: (1, 1), 2: (1, 0), 3: (1, -1)}

CLOUD_DTYPE = np.dtype([
    ("family", "u1"), ("row", "f8"), ("col", "f8"), ("depth_m", "f8"), ("weight", "f8"),
    ("depth_se_m", "f8"), ("residual", "f8"),
])


def load_cloud_csv(path: str | Path, transform) -> np.ndarray:
    """Load the frozen H4 solution cloud; keep only cluster-weighted records.

    ``transform`` must be the exact sample submission affine (100 m, north-up).
    """
    a, _, c, _, e, f = tuple(transform)[:6]
    cloud = []
    with gzip.open(path, "rb") as raw:
        with io.TextIOWrapper(raw, encoding="utf-8", newline="") as text:
            reader = csv.reader(text)
            header = next(reader)
            required = ["family", "easting_m", "northing_m", "effective_depth_m",
                        "conditional_depth_se_m", "relative_residual", "cluster_weight"]
            index = {name: header.index(name) for name in required}
            n_rows = 0
            for row in reader:
                n_rows += 1
                weight = float(row[index["cluster_weight"]])
                if weight <= 0.0:
                    continue
                east = float(row[index["easting_m"]])
                north = float(row[index["northing_m"]])
                col = (east - c) / a - 0.5
                rowf = (north - f) / e - 0.5
                family = FAMILIES.get(row[index["family"]])
                if family is None or not all(np.isfinite([rowf, col])):
                    continue
                cloud.append((family, rowf, col, float(row[index["effective_depth_m"]]), weight,
                              float(row[index["conditional_depth_se_m"]]),
                              float(row[index["relative_residual"]])))
    if n_rows != N_RECORDS:
        raise ValueError(f"cloud record count {n_rows} != pinned {N_RECORDS}")
    if not cloud:
        raise ValueError("no cluster-weighted solutions in the frozen cloud")
    return np.array(cloud, dtype=CLOUD_DTYPE)


def _kde(rows: np.ndarray, cols: np.ndarray, weights: np.ndarray, footprint: np.ndarray) -> np.ndarray:
    raw = bilinear_splat(rows, cols, weights, footprint.shape)
    kde = ndimage.gaussian_filter(raw, sigma=KDE_SIGMA_PX, truncate=KDE_TRUNCATE,
                                  mode="constant", cval=0.0)
    kde[~footprint] = 0.0
    support = kde[footprint & (kde > 0)]
    if not len(support):
        raise ValueError("no KDE support in footprint; stop without relaxing gates")
    scale = float(np.percentile(support, FAMILY_PERCENTILE))
    return np.clip(kde / scale, 0.0, 1.0).astype(np.float32)


def family_kdes(cloud: np.ndarray, footprint: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict]:
    mag = cloud[cloud["family"] == FAMILIES["tmi"]]
    grav = cloud[cloud["family"] == FAMILIES["iso_grav_anom"]]
    if not len(mag) or not len(grav):
        raise ValueError("both potential-field families must contribute weighted solutions")
    m = _kde(mag["row"], mag["col"], mag["weight"], footprint)
    g = _kde(grav["row"], grav["col"], grav["weight"], footprint)
    stats = {"magnetic_weighted_solutions": int(len(mag)), "gravity_weighted_solutions": int(len(grav)),
             "magnetic_positive_cells": int(np.count_nonzero(m)), "gravity_positive_cells": int(np.count_nonzero(g)),
             "magnetic_depth_p50_m": float(np.median(mag["depth_m"])),
             "gravity_depth_p50_m": float(np.median(grav["depth_m"]))}
    return m, g, stats


def consensus_field(magnetic: np.ndarray, gravity: np.ndarray, footprint: np.ndarray,
                    known: np.ndarray) -> np.ndarray:
    """Frozen H4 combine rule, with exact known pixels zeroed before orientation analysis."""
    field = (magnetic + gravity + np.sqrt(magnetic * gravity)) / 3.0
    field[~footprint | known] = 0.0
    return field.astype(np.float32)


def orientation_bins(consensus: np.ndarray, footprint: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Quantized gradient direction (mod 180) and anisotropy of the consensus field.

    Returns ``(bins, aniso)`` with bins in {0..3} where oriented and -1 elsewhere.
    """
    inner = ndimage.gaussian_filter(consensus, sigma=TENSOR_SIGMA_INNER_PX, mode="constant", cval=0.0)
    gy, gx = np.gradient(inner)
    integration = TENSOR_SIGMA_INTEGRATION_PX
    jxx = ndimage.gaussian_filter(gx * gx, sigma=integration, mode="constant", cval=0.0)
    jyy = ndimage.gaussian_filter(gy * gy, sigma=integration, mode="constant", cval=0.0)
    jxy = ndimage.gaussian_filter(gx * gy, sigma=integration, mode="constant", cval=0.0)
    trace = jxx + jyy
    det = jxx * jyy - jxy * jxy
    disc = np.sqrt(np.maximum(trace * trace / 4.0 - det, 0.0))
    lam1, lam2 = trace / 2.0 + disc, trace / 2.0 - disc
    aniso = np.where(trace > 0, (lam1 - lam2) / np.maximum(trace, 1e-30), 0.0)
    alpha = 0.5 * np.arctan2(2.0 * jxy, jxx - jyy)  # gradient direction from the column axis
    bins = np.round(alpha / (np.pi / 4.0)).astype(np.int8) % 4
    oriented = footprint & (consensus > 0) & (aniso >= ISOTROPY_MIN)
    bins = np.where(oriented, bins, -1).astype(np.int8)
    return bins, aniso.astype(np.float32)


def _shift(arr: np.ndarray, dr: int, dc: int) -> np.ndarray:
    """``out[r, c] = arr[r + dr, c + dc]`` with zero fill (no wrap-around)."""
    height, width = arr.shape
    out = np.zeros_like(arr)
    r0, r1 = max(0, -dr), min(height, height - dr)
    c0, c1 = max(0, -dc), min(width, width - dc)
    out[r0:r1, c0:c1] = arr[r0 + dr:r1 + dr, c0 + dc:c1 + dc]
    return out


def ridge_crests(consensus: np.ndarray, bins: np.ndarray) -> np.ndarray:
    """Orientation-NMS crests: local maxima along the quantized gradient direction.

    Plateau-safe one-sided tie rule: ``C >= forward`` and ``C > backward`` thins constant
    plateaus to a single 1-px ridge regardless of plateau width.
    """
    ridge = np.zeros(consensus.shape, dtype=bool)
    for b, (dr, dc) in BIN_DIRECTIONS.items():
        mask = bins == b
        if not mask.any():
            continue
        fwd = _shift(consensus, dr, dc)    # neighbour one step along +gradient
        bwd = _shift(consensus, -dr, -dc)  # neighbour one step against +gradient
        ridge |= mask & (consensus > 0) & (consensus >= fwd) & (consensus > bwd)
    return ridge


def corridor_kernels() -> dict[int, np.ndarray]:
    """Oriented along-strike box kernels, one per gradient bin (strike = gradient rotated 90)."""
    half_l, half_w = CORRIDOR_ALONG_HALF_PX, CORRIDOR_ACROSS_HALF_PX
    size = 2 * half_l + 1
    kernels = {}
    for b, (dr, dc) in BIN_DIRECTIONS.items():
        norm = float(np.hypot(dr, dc))
        grd = np.array((dr, dc), dtype=float) / norm
        strike = np.array((-dc, dr), dtype=float)  # perpendicular
        kernel = np.zeros((size, size), dtype=np.float32)
        centre = half_l
        for i in range(size):
            for j in range(size):
                vec = np.array((i - centre, j - centre), dtype=float)
                along = float(np.dot(vec, strike))
                across = float(np.dot(vec, grd))
                if abs(along) <= half_l + 0.5 and abs(across) <= half_w + 0.5:
                    kernel[i, j] = 1.0
        if not kernel.any():
            raise ValueError("empty corridor kernel")
        kernels[b] = kernel
    return kernels


def strike_consensus(cloud: np.ndarray, bins: np.ndarray, ridge: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict]:
    """Corridor solution count and depth dispersion per crest cell.

    Returns ``(passes, gamma, stats)``.  ``passes`` marks crests whose corridor contains
    >= MIN_CORRIDOR_SOLUTIONS solutions (raw count) whose *cluster-weight-weighted* depth
    standard deviation <= CORRIDOR_SD_MAX_M.  The moments are weight-weighted because the frozen
    H4 cluster weight already down-weights deep, scattered and mutually inconsistent solutions;
    counting those equally would let incoherent noise veto coherent shallow clusters
    (pre-scoring implementation correction recorded in the preregistration, section 6).
    """
    shape = ridge.shape
    w = cloud["weight"]
    count_grid = bilinear_splat(cloud["row"], cloud["col"], np.ones_like(w), shape)
    wsum_grid = bilinear_splat(cloud["row"], cloud["col"], w, shape)
    wz_grid = bilinear_splat(cloud["row"], cloud["col"], w * cloud["depth_m"], shape)
    wz2_grid = bilinear_splat(cloud["row"], cloud["col"], w * cloud["depth_m"] ** 2, shape)
    n = np.zeros(shape, dtype=np.float64)
    mu = np.zeros(shape, dtype=np.float64)
    sd = np.full(shape, np.inf, dtype=np.float64)
    for b, kernel in corridor_kernels().items():
        sel = (bins == b) & ridge
        if not sel.any():
            continue
        nb = ndimage.convolve(count_grid, kernel, mode="constant", cval=0.0)
        wb = ndimage.convolve(wsum_grid, kernel, mode="constant", cval=0.0)
        wzb = ndimage.convolve(wz_grid, kernel, mode="constant", cval=0.0)
        wz2b = ndimage.convolve(wz2_grid, kernel, mode="constant", cval=0.0)
        with np.errstate(invalid="ignore", divide="ignore"):
            mu_b = np.where(wb > 0, wzb / wb, 0.0)
            var_b = np.where(wb > 0, np.maximum(wz2b / wb - mu_b ** 2, 0.0), np.inf)
        n[sel] = nb[sel]
        mu[sel] = mu_b[sel]
        sd[sel] = np.sqrt(var_b)[sel]
    passes = ridge & (n >= MIN_CORRIDOR_SOLUTIONS) & (sd <= CORRIDOR_SD_MAX_M)
    gamma = np.where(passes, np.exp(-np.where(np.isfinite(sd), sd, CORRIDOR_SD_MAX_M) / CORRIDOR_SD_MAX_M), 0.0)
    stats = {"ridge_cells": int(ridge.sum()), "corridor_passing_cells": int(passes.sum()),
             "passing_fraction": float(passes.sum() / ridge.sum()) if ridge.any() else 0.0,
             "passing_median_depth_m": float(np.median(mu[passes])) if passes.any() else None,
             "passing_median_sd_m": float(np.median(sd[passes])) if passes.any() else None}
    return passes, gamma.astype(np.float32), stats


def build_trace_locked_field(cloud: np.ndarray, footprint: np.ndarray, known: np.ndarray) -> tuple[np.ndarray, dict]:
    """Full H8 construction; fails closed rather than relaxing any frozen gate."""
    m, g, kde_stats = family_kdes(cloud, footprint)
    consensus = consensus_field(m, g, footprint, known)
    bins, aniso = orientation_bins(consensus, footprint)
    ridge = ridge_crests(consensus, bins)
    passes, gamma, stats = strike_consensus(cloud, bins, ridge)
    if not passes.any():
        raise ValueError("no crest passed the along-strike depth consensus; stop, do not relax gates")
    raw = np.where(passes, consensus * gamma, 0.0).astype(np.float32)
    field = ndimage.gaussian_filter(raw, sigma=OUTPUT_SIGMA_PX, truncate=OUTPUT_TRUNCATE,
                                    mode="constant", cval=0.0)
    field[~footprint] = 0.0
    maximum = float(field[footprint].max())
    if not np.isfinite(maximum) or maximum <= 0.0:
        raise ValueError("empty trace-locked field")
    field = np.clip(field / maximum, 0.0, 1.0).astype(np.float32)
    field[known] = 0.0
    stats = {**kde_stats, **stats,
             "oriented_cells": int((bins >= 0).sum()),
             "positive_output_cells": int(np.count_nonzero(field[footprint])),
             "known_pixels_zeroed": int(known.sum()),
             "output_mass_sum": float(field[footprint].sum())}
    return field, stats
