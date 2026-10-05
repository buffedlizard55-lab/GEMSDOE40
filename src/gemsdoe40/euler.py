"""Preregistered H1: paired magnetic/gravity 3-D Euler clouds and depth-consensus KDE."""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy import ndimage
from scipy import fft as sp_fft
from scipy.spatial import cKDTree
import rasterio

from .raster import band_index_by_name


WINDOW = 10
STRIDE = 5
RESOLUTION_M = 100.0
MAX_RESIDUAL = 0.20
MAX_CONDITION = 1e5
MAX_HORIZONTAL_OFFSET_M = WINDOW * RESOLUTION_M / 2.0
MIN_DEPTH_M = 100.0
MAX_DEPTH_M = 10_000.0
PAIR_MAX_XY_M = 300.0
PAIR_MAX_DEPTH_M = 1_000.0
KDE_SIGMA_PX = 2.0
KDE_TRUNCATE = 3.0
EMISSION_BUDGET = 45_962
FOURIER_PAD_CELLS = 128


@dataclass
class EulerCloud:
    row: np.ndarray
    col: np.ndarray
    depth_m: np.ndarray
    residual: np.ndarray
    condition: np.ndarray

    def __len__(self) -> int:
        return int(self.row.size)

    def to_summary(self) -> dict[str, Any]:
        if not len(self):
            return {"n_solutions": 0}
        return {
            "n_solutions": len(self),
            "depth_m_p05_p50_p95": np.percentile(self.depth_m, [5, 50, 95]).tolist(),
            "residual_p50_p95": np.percentile(self.residual, [50, 95]).tolist(),
            "condition_p50_p95": np.percentile(self.condition, [50, 95]).tolist(),
        }


@dataclass
class PairCloud:
    row: np.ndarray
    col: np.ndarray
    weight: np.ndarray
    xy_distance_m: np.ndarray
    depth_delta_m: np.ndarray

    def to_summary(self) -> dict[str, Any]:
        if not self.row.size:
            return {"n_pairs": 0}
        return {
            "n_pairs": int(self.row.size),
            "pair_weight_sum": float(self.weight.sum(dtype=np.float64)),
            "xy_distance_m_p50_p95": np.percentile(self.xy_distance_m, [50, 95]).tolist(),
            "depth_delta_m_p50_p95": np.percentile(self.depth_delta_m, [50, 95]).tolist(),
        }


def _finite_mask(array: np.ndarray) -> np.ndarray:
    arr = np.asarray(array)
    return np.isfinite(arr) & (np.abs(arr) < 1e30)


def fill_nearest(array: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Fill nodata only for derivative calculations; return the original valid mask."""
    arr = np.asarray(array, dtype=np.float32)
    valid = _finite_mask(arr)
    if not valid.any():
        raise ValueError("cannot interpolate a band with no finite values")
    if valid.all():
        return arr.copy(), valid
    nearest = ndimage.distance_transform_edt(~valid, return_distances=False, return_indices=True)
    filled = arr.copy()
    filled[~valid] = arr[tuple(nearest)][~valid]
    return filled, valid


def fourier_vertical_down(array: np.ndarray, *, pad_cells: int = FOURIER_PAD_CELLS, resolution_m: float = RESOLUTION_M) -> np.ndarray:
    """Compute a downward-positive vertical derivative (+|k|) of a harmonic field.

    A reflect-padded grid reduces rectangular-edge ringing. Invalid cells are
    nearest-filled strictly for this derivative calculation and are excluded
    from Euler windows by the caller.
    """
    filled, _ = fill_nearest(array)
    padded = np.pad(filled, ((pad_cells, pad_cells), (pad_cells, pad_cells)), mode="reflect")
    ny, nx = padded.shape
    fy = sp_fft.fftfreq(ny, d=resolution_m).astype(np.float32)[:, None]
    fx = sp_fft.rfftfreq(nx, d=resolution_m).astype(np.float32)[None, :]
    wave_number = (2.0 * np.pi * np.sqrt(fy * fy + fx * fx)).astype(np.float32)
    spectrum = sp_fft.rfft2(padded, workers=-1)
    deriv = sp_fft.irfft2(spectrum * wave_number, s=padded.shape, workers=-1)
    deriv = deriv[pad_cells:-pad_cells, pad_cells:-pad_cells]
    return np.asarray(deriv, dtype=np.float32)


def upward_continue_with_vertical(
    array: np.ndarray,
    height_m: float,
    *,
    pad_cells: int = FOURIER_PAD_CELLS,
    resolution_m: float = RESOLUTION_M,
) -> tuple[np.ndarray, np.ndarray]:
    """Upward-continue a harmonic field and derive its downward-positive Tz."""
    if height_m < 0:
        raise ValueError("upward-continuation height must be non-negative")
    filled, _ = fill_nearest(array)
    padded = np.pad(filled, ((pad_cells, pad_cells), (pad_cells, pad_cells)), mode="reflect")
    ny, nx = padded.shape
    fy = sp_fft.fftfreq(ny, d=resolution_m).astype(np.float32)[:, None]
    fx = sp_fft.rfftfreq(nx, d=resolution_m).astype(np.float32)[None, :]
    wave_number = (2.0 * np.pi * np.sqrt(fy * fy + fx * fx)).astype(np.float32)
    spectrum = sp_fft.rfft2(padded, workers=-1)
    continue_filter = np.exp(-wave_number * np.float32(height_m)).astype(np.float32)
    continued_spectrum = spectrum * continue_filter
    continued = sp_fft.irfft2(continued_spectrum, s=padded.shape, workers=-1)
    vertical = sp_fft.irfft2(continued_spectrum * wave_number, s=padded.shape, workers=-1)
    crop = (slice(pad_cells, -pad_cells), slice(pad_cells, -pad_cells))
    return np.asarray(continued[crop], dtype=np.float32), np.asarray(vertical[crop], dtype=np.float32)


def _gradient_components(array: np.ndarray, resolution_m: float = RESOLUTION_M) -> tuple[np.ndarray, np.ndarray]:
    """Return derivatives in positive-east and positive-north coordinates."""
    dy_row, dx_col = np.gradient(np.asarray(array, dtype=np.float32), resolution_m, edge_order=1)
    return dx_col.astype(np.float32, copy=False), (-dy_row).astype(np.float32, copy=False)


def validate_magnetic_vertical_derivative(tmi: np.ndarray, supplied_tmi_vg: np.ndarray) -> dict[str, float | bool]:
    """Compare the named TMI vertical-gradient feature with +|k| derived TMI."""
    spectral = fourier_vertical_down(tmi)
    supplied = np.asarray(supplied_tmi_vg, dtype=np.float32)
    good = _finite_mask(tmi) & _finite_mask(supplied)
    # Avoid the outermost stencil rows/columns; the sample mask is applied later.
    good[:5, :] = False
    good[-5:, :] = False
    good[:, :5] = False
    good[:, -5:] = False
    x = spectral[good].astype(np.float64)
    y = supplied[good].astype(np.float64)
    if x.size < 100:
        raise ValueError("not enough finite cells to validate tmi_vg")
    corr = float(np.corrcoef(x, y)[0, 1])
    denom = float(np.dot(x, x))
    scale = float(np.dot(x, y) / denom) if denom > 0 else float("nan")
    return {"pearson_r": corr, "supplied_over_spectral_scale": scale, "pass": bool(corr >= 0.95 and 0.8 <= scale <= 1.2)}


def solve_euler_patch(
    dtx: np.ndarray,
    dty: np.ndarray,
    dtz: np.ndarray,
    *,
    resolution_m: float = RESOLUTION_M,
) -> tuple[float, float, float, float, float]:
    """Solve one SI=0 patch, returning x0/y0 offsets, depth, residual, condition.

    Coordinates are centered on the window in metres; y is positive north and
    depth is positive down. For SI=0, Euler homogeneity gives
    A @ [x0, y0, z0] = x*Tx + y*Ty, since the observation plane has z=0.
    """
    dtx = np.asarray(dtx, dtype=np.float64)
    dty = np.asarray(dty, dtype=np.float64)
    dtz = np.asarray(dtz, dtype=np.float64)
    if dtx.shape != (WINDOW, WINDOW) or dty.shape != dtx.shape or dtz.shape != dtx.shape:
        raise ValueError(f"Euler patch must be {WINDOW}x{WINDOW}")
    if not (np.isfinite(dtx).all() and np.isfinite(dty).all() and np.isfinite(dtz).all()):
        raise ValueError("Euler patch derivatives must be finite")
    offsets = (np.arange(WINDOW, dtype=np.float64) + 0.5 - WINDOW / 2.0) * resolution_m
    x_obs = np.broadcast_to(offsets[None, :], dtx.shape)
    y_obs = np.broadcast_to((-offsets)[:, None], dtx.shape)
    design = np.stack((dtx.ravel(), dty.ravel(), dtz.ravel()), axis=1)
    rhs = (x_obs * dtx + y_obs * dty).ravel()
    normal = design.T @ design
    cond = float(np.linalg.cond(normal))
    if not np.isfinite(cond) or cond > MAX_CONDITION:
        return (float("nan"),) * 5
    cross = design.T @ rhs
    solution = np.linalg.solve(normal, cross)
    residual_sq = float(np.dot(rhs, rhs) - np.dot(cross, solution))
    relative_residual = math.sqrt(max(0.0, residual_sq) / (float(np.dot(rhs, rhs)) + 1e-30))
    return float(solution[0]), float(solution[1]), float(solution[2]), relative_residual, cond


def solve_euler_cloud(
    tx: np.ndarray,
    ty: np.ndarray,
    tz: np.ndarray,
    window_valid: np.ndarray,
    *,
    stride: int = STRIDE,
    chunk_window_rows: int = 8,
    window: int = WINDOW,
    resolution_m: float = RESOLUTION_M,
    height_m: float = 0.0,
) -> EulerCloud:
    """Solve sliding, overlapping SI=0 Euler windows and apply locked QC gates."""
    tx, ty, tz = (np.asarray(x, dtype=np.float32) for x in (tx, ty, tz))
    valid = np.asarray(window_valid, dtype=bool)
    if tx.ndim != 2 or tx.shape != ty.shape or tx.shape != tz.shape or tx.shape != valid.shape:
        raise ValueError("derivative arrays and window_valid must have one identical 2-D shape")
    if window != WINDOW or stride <= 0 or chunk_window_rows <= 0:
        raise ValueError("H1 uses a 10-cell window; stride and chunk size must be positive")
    h, w = tx.shape
    if h < window or w < window:
        return EulerCloud(*(np.empty(0, dtype=np.float32) for _ in range(5)))

    # Require all 10x10 cells and the one-cell central-difference halo to have
    # real input support. The 3x3 minimum filter applies the halo requirement.
    deriv_valid = ndimage.minimum_filter(valid.astype(np.uint8), size=3, mode="constant", cval=0) > 0
    valid_windows = np.lib.stride_tricks.sliding_window_view(deriv_valid, (window, window))[::stride, ::stride]
    valid_windows = valid_windows.all(axis=(-1, -2))
    txw = np.lib.stride_tricks.sliding_window_view(tx, (window, window))[::stride, ::stride]
    tyw = np.lib.stride_tricks.sliding_window_view(ty, (window, window))[::stride, ::stride]
    tzw = np.lib.stride_tricks.sliding_window_view(tz, (window, window))[::stride, ::stride]
    nyw, nxw = valid_windows.shape
    xoff = (np.arange(window, dtype=np.float32) + 0.5 - window / 2.0) * resolution_m
    xgrid = np.broadcast_to(xoff[None, :], (window, window))
    ygrid = np.broadcast_to((-xoff)[:, None], (window, window))

    rows_out: list[np.ndarray] = []
    cols_out: list[np.ndarray] = []
    depths_out: list[np.ndarray] = []
    residuals_out: list[np.ndarray] = []
    conditions_out: list[np.ndarray] = []
    for r0 in range(0, nyw, chunk_window_rows):
        r1 = min(nyw, r0 + chunk_window_rows)
        # Each chunk stays bounded in memory; these are strided read-only views.
        ax, ay, az = txw[r0:r1], tyw[r0:r1], tzw[r0:r1]
        design = np.stack((ax, ay, az), axis=-1).reshape(r1 - r0, nxw, window * window, 3).astype(np.float64)
        rhs_image = xgrid * ax + ygrid * ay
        rhs = rhs_image.reshape(r1 - r0, nxw, window * window).astype(np.float64)
        normal = np.einsum("...ki,...kj->...ij", design, design, optimize=True)
        cross = np.einsum("...ki,...k->...i", design, rhs, optimize=True)
        rhs_norm = np.einsum("...k,...k->...", rhs, rhs, optimize=True)
        eigenvalues = np.linalg.eigvalsh(normal)
        cond = np.full(eigenvalues.shape[:-1], np.inf, dtype=np.float64)
        eig_ok = eigenvalues[..., 0] > 0
        cond[eig_ok] = eigenvalues[..., -1][eig_ok] / eigenvalues[..., 0][eig_ok]
        good = valid_windows[r0:r1] & np.isfinite(cond) & (cond <= MAX_CONDITION) & (rhs_norm > 1e-24)
        if not good.any():
            continue
        # Solve only well-conditioned patches. A batched solve is much faster
        # than a Python loop but still uses normal equations only after QC.
        normal_good = normal[good]
        cross_good = cross[good]
        try:
            solutions = np.linalg.solve(normal_good, cross_good[..., None])[..., 0]
        except np.linalg.LinAlgError:
            # Defensive per-window fallback for rare numerical singularities.
            solved: list[np.ndarray] = []
            good_locations = np.argwhere(good)
            for (rr, cc), mat, vec in zip(good_locations, normal_good, cross_good):
                try:
                    solved.append(np.linalg.solve(mat, vec))
                except np.linalg.LinAlgError:
                    solved.append(np.full(3, np.nan))
            solutions = np.asarray(solved, dtype=np.float64)
        fit = np.einsum("...i,...i->...", cross_good, solutions, optimize=True)
        sse = np.maximum(rhs_norm[good] - fit, 0.0)
        residual = np.sqrt(sse / (rhs_norm[good] + 1e-30))
        loc = np.argwhere(good)
        start_rows = (loc[:, 0] + r0) * stride
        start_cols = loc[:, 1] * stride
        dx0, dy0, depth_from_observation = solutions.T
        depth_ground = depth_from_observation - height_m
        source_col = start_cols + (window / 2.0 - 0.5) + dx0 / resolution_m
        source_row = start_rows + (window / 2.0 - 0.5) - dy0 / resolution_m
        cond_good = cond[good]
        accept = (
            np.isfinite(source_col)
            & np.isfinite(source_row)
            & np.isfinite(depth_ground)
            & (residual <= MAX_RESIDUAL)
            & (np.abs(dx0) <= MAX_HORIZONTAL_OFFSET_M)
            & (np.abs(dy0) <= MAX_HORIZONTAL_OFFSET_M)
            & (depth_ground >= MIN_DEPTH_M)
            & (depth_ground <= MAX_DEPTH_M)
        )
        if accept.any():
            rows_out.append(source_row[accept].astype(np.float32))
            cols_out.append(source_col[accept].astype(np.float32))
            depths_out.append(depth_ground[accept].astype(np.float32))
            residuals_out.append(residual[accept].astype(np.float32))
            conditions_out.append(cond_good[accept].astype(np.float32))

    def combine(parts: list[np.ndarray]) -> np.ndarray:
        return np.concatenate(parts) if parts else np.empty(0, dtype=np.float32)

    return EulerCloud(*(combine(parts) for parts in (rows_out, cols_out, depths_out, residuals_out, conditions_out)))


def pair_clouds(magnetic: EulerCloud, gravity: EulerCloud) -> PairCloud:
    """Pair each magnetic solution to its best nearby depth-consistent gravity solution."""
    if not len(magnetic) or not len(gravity):
        empty = np.empty(0, dtype=np.float32)
        return PairCloud(empty, empty, empty, empty, empty)
    grav_xy = np.column_stack((gravity.col * RESOLUTION_M, gravity.row * RESOLUTION_M)).astype(np.float64)
    mag_xy = np.column_stack((magnetic.col * RESOLUTION_M, magnetic.row * RESOLUTION_M)).astype(np.float64)
    tree = cKDTree(grav_xy)
    pair_rows: list[float] = []
    pair_cols: list[float] = []
    pair_weights: list[float] = []
    pair_distances: list[float] = []
    pair_depth_deltas: list[float] = []
    chunk = 4096
    for start in range(0, len(magnetic), chunk):
        stop = min(len(magnetic), start + chunk)
        neighbors = tree.query_ball_point(mag_xy[start:stop], r=PAIR_MAX_XY_M, workers=-1, return_sorted=True)
        for local, candidates in enumerate(neighbors):
            if not candidates:
                continue
            mi = start + local
            gi = np.asarray(candidates, dtype=np.int64)
            dxy = np.linalg.norm(grav_xy[gi] - mag_xy[mi], axis=1)
            dz = np.abs(gravity.depth_m[gi] - magnetic.depth_m[mi])
            valid = dz <= PAIR_MAX_DEPTH_M
            if not valid.any():
                continue
            gi, dxy, dz = gi[valid], dxy[valid], dz[valid]
            joint = (dxy / PAIR_MAX_XY_M) ** 2 + (dz / PAIR_MAX_DEPTH_M) ** 2
            best = int(np.argmin(joint))
            gj = int(gi[best])
            distance, depth_delta = float(dxy[best]), float(dz[best])
            residual_term = (float(magnetic.residual[mi]) / MAX_RESIDUAL) ** 2 + (float(gravity.residual[gj]) / MAX_RESIDUAL) ** 2
            spatial_term = (distance / PAIR_MAX_XY_M) ** 2 + (depth_delta / PAIR_MAX_DEPTH_M) ** 2
            weight = math.exp(-0.5 * (residual_term + spatial_term))
            pair_rows.append(0.5 * (float(magnetic.row[mi]) + float(gravity.row[gj])))
            pair_cols.append(0.5 * (float(magnetic.col[mi]) + float(gravity.col[gj])))
            pair_weights.append(weight)
            pair_distances.append(distance)
            pair_depth_deltas.append(depth_delta)
    return PairCloud(
        row=np.asarray(pair_rows, dtype=np.float32),
        col=np.asarray(pair_cols, dtype=np.float32),
        weight=np.asarray(pair_weights, dtype=np.float32),
        xy_distance_m=np.asarray(pair_distances, dtype=np.float32),
        depth_delta_m=np.asarray(pair_depth_deltas, dtype=np.float32),
    )


def persistence_cloud(clouds_by_height: dict[int, EulerCloud]) -> tuple[PairCloud, np.ndarray]:
    """Retain height-0 solutions that persist across at least 3/4 heights."""
    heights = [0, 500, 1000, 2000]
    if any(height not in clouds_by_height for height in heights):
        raise ValueError(f"persistence clouds must include heights {heights}")
    base = clouds_by_height[0]
    if not len(base):
        empty = np.empty(0, dtype=np.float32)
        return PairCloud(empty, empty, empty, empty, empty), empty.astype(np.int8)
    trees: dict[int, cKDTree] = {}
    xy_by_height: dict[int, np.ndarray] = {}
    for height in heights[1:]:
        cloud = clouds_by_height[height]
        xy = np.column_stack((cloud.col * RESOLUTION_M, cloud.row * RESOLUTION_M)).astype(np.float64)
        xy_by_height[height] = xy
        if len(cloud):
            trees[height] = cKDTree(xy)

    rows_out: list[float] = []
    cols_out: list[float] = []
    weights_out: list[float] = []
    mean_xy_out: list[float] = []
    mean_depth_out: list[float] = []
    hits_out: list[int] = []
    base_xy = np.column_stack((base.col * RESOLUTION_M, base.row * RESOLUTION_M)).astype(np.float64)
    for mi in range(len(base)):
        hits = 1
        terms = [(float(base.residual[mi]) / MAX_RESIDUAL) ** 2]
        xy_values: list[float] = []
        depth_values: list[float] = []
        for height in heights[1:]:
            cloud = clouds_by_height[height]
            tree = trees.get(height)
            if tree is None:
                continue
            candidates = tree.query_ball_point(base_xy[mi], r=PAIR_MAX_XY_M, workers=-1, return_sorted=True)
            if not candidates:
                continue
            gi = np.asarray(candidates, dtype=np.int64)
            dxy = np.linalg.norm(xy_by_height[height][gi] - base_xy[mi], axis=1)
            dz = np.abs(cloud.depth_m[gi] - base.depth_m[mi])
            allowed = dz <= PAIR_MAX_DEPTH_M
            if not allowed.any():
                continue
            gi, dxy, dz = gi[allowed], dxy[allowed], dz[allowed]
            joint = (dxy / PAIR_MAX_XY_M) ** 2 + (dz / PAIR_MAX_DEPTH_M) ** 2
            best = int(np.argmin(joint))
            gj = int(gi[best])
            hits += 1
            terms.append(
                (float(clouds_by_height[height].residual[gj]) / MAX_RESIDUAL) ** 2
                + float(joint[best])
            )
            xy_values.append(float(dxy[best]))
            depth_values.append(float(dz[best]))
        if hits < 3:
            continue
        mean_term = float(np.mean(terms))
        weight = (hits / 4.0) * math.exp(-0.5 * mean_term)
        rows_out.append(float(base.row[mi]))
        cols_out.append(float(base.col[mi]))
        weights_out.append(weight)
        mean_xy_out.append(float(np.mean(xy_values)) if xy_values else 0.0)
        mean_depth_out.append(float(np.mean(depth_values)) if depth_values else 0.0)
        hits_out.append(hits)
    cloud = PairCloud(
        row=np.asarray(rows_out, dtype=np.float32),
        col=np.asarray(cols_out, dtype=np.float32),
        weight=np.asarray(weights_out, dtype=np.float32),
        xy_distance_m=np.asarray(mean_xy_out, dtype=np.float32),
        depth_delta_m=np.asarray(mean_depth_out, dtype=np.float32),
    )
    return cloud, np.asarray(hits_out, dtype=np.int8)


def kde_from_pairs(pairs: PairCloud, shape: tuple[int, int]) -> np.ndarray:
    """Deposit pair weights, then apply the preregistered zero-padded Gaussian KDE."""
    density = np.zeros(shape, dtype=np.float32)
    if pairs.row.size:
        rr = np.rint(pairs.row).astype(np.int64)
        cc = np.rint(pairs.col).astype(np.int64)
        inside = (rr >= 0) & (rr < shape[0]) & (cc >= 0) & (cc < shape[1])
        np.add.at(density, (rr[inside], cc[inside]), pairs.weight[inside])
    return ndimage.gaussian_filter(density, sigma=KDE_SIGMA_PX, truncate=KDE_TRUNCATE, mode="constant", cval=0.0)


def emit_top_budget(
    kde: np.ndarray,
    valid: np.ndarray,
    known_labels: np.ndarray,
    *,
    budget: int = EMISSION_BUDGET,
) -> tuple[np.ndarray, dict[str, int | float]]:
    """Emit binary top-budget cells with deterministic row-major tie breaking."""
    kde = np.asarray(kde, dtype=np.float32)
    valid = np.asarray(valid, dtype=bool)
    known = np.asarray(known_labels)
    if kde.shape != valid.shape or known.shape != kde.shape:
        raise ValueError("KDE, valid mask, and known labels must have the same shape")
    eligible = valid & (known != 1) & np.isfinite(kde) & (kde > 0.0)
    indices = np.flatnonzero(eligible)
    if indices.size < budget:
        raise RuntimeError(f"only {indices.size} positive KDE cells; preregistered budget is {budget}; refuse to backfill")
    values = kde.ravel()[indices]
    kth = float(np.partition(values, values.size - budget)[values.size - budget])
    higher = indices[values > kth]
    equal = np.sort(indices[values == kth])
    need = budget - higher.size
    selected = np.concatenate((higher, equal[:need]))
    if selected.size != budget:
        raise AssertionError(f"selected {selected.size} cells, expected {budget}")
    unselected_ties = equal[need:]
    if unselected_ties.size:
        maximum_unemitted = kth
    else:
        below = values[values < kth]
        maximum_unemitted = float(below.max()) if below.size else 0.0
    pred = np.zeros(kde.shape, dtype=np.float32)
    pred.ravel()[selected] = 1.0
    summary = {
        "budget": int(budget),
        "positive_kde_cells": int(indices.size),
        "emitted_pixels": int(selected.size),
        "cutoff_kde": kth,
        "minimum_emitted_kde": float(kde.ravel()[selected].min()),
        "maximum_unemitted_kde": maximum_unemitted,
    }
    return pred, summary


def emit_all_supported(
    kde: np.ndarray,
    valid: np.ndarray,
    known_labels: np.ndarray,
) -> tuple[np.ndarray, dict[str, int | float | str]]:
    """H2-B: emit exactly the finite positive KDE support; no budget fill/tuning."""
    kde = np.asarray(kde, dtype=np.float32)
    valid = np.asarray(valid, dtype=bool)
    known = np.asarray(known_labels)
    if kde.shape != valid.shape or known.shape != kde.shape:
        raise ValueError("KDE, valid mask, and known labels must have the same shape")
    eligible = valid & (known != 1) & np.isfinite(kde) & (kde > 0.0)
    count = int(eligible.sum())
    if count == 0:
        raise RuntimeError("H2-B has no positive KDE support; no candidate can be emitted")
    pred = np.zeros(kde.shape, dtype=np.float32)
    pred[eligible] = 1.0
    return pred, {
        "emission_rule": "all positive KDE support; no top-K threshold or backfill",
        "positive_kde_cells": count,
        "emitted_pixels": count,
        "minimum_emitted_kde": float(kde[eligible].min()),
        "maximum_unemitted_kde": 0.0,
    }


def run_h1(
    features_path: str | Path,
    template_path: str | Path,
    labels_path: str | Path,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Generate one frozen H1 candidate from challenge feature bands only."""
    features_path, template_path, labels_path = Path(features_path), Path(template_path), Path(labels_path)
    with rasterio.open(template_path) as sample, rasterio.open(labels_path) as labels_ds, rasterio.open(features_path) as features:
        template = sample.read(1)
        footprint = np.isfinite(template)
        labels = labels_ds.read(1)
        if labels.shape != template.shape or features.shape != template.shape:
            raise ValueError("sample, labels, and features must have identical shape")
        bands: dict[str, np.ndarray] = {}
        band_indices: dict[str, int] = {}
        for name in ("tmi", "tmi_vg", "iso_grav_anom"):
            idx = band_index_by_name(features, name)
            band_indices[name] = idx
            bands[name] = features.read(idx).astype(np.float32, copy=False)
        transform = features.transform
        if sample.crs != features.crs or sample.transform != features.transform:
            raise ValueError("sample and feature georeferencing differ")

    tmi, tmi_mask = fill_nearest(bands["tmi"])
    tmi_vg, tmi_vg_mask = fill_nearest(bands["tmi_vg"])
    grav, grav_mask = fill_nearest(bands["iso_grav_anom"])
    common_valid = footprint & tmi_mask & tmi_vg_mask & grav_mask
    derivative_check = validate_magnetic_vertical_derivative(bands["tmi"], bands["tmi_vg"])
    if derivative_check["pass"]:
        magnetic_tz = tmi_vg
        magnetic_vertical_source = "provided tmi_vg, independently verified against +|k| TMI"
    else:
        magnetic_tz = fourier_vertical_down(tmi)
        magnetic_vertical_source = "Fourier +|k| derivative fallback (tmi_vg failed preregistered check)"

    magnetic_tx, magnetic_ty = _gradient_components(tmi)
    gravity_tz = fourier_vertical_down(grav)
    gravity_tx, gravity_ty = _gradient_components(gravity_tz)
    gravity_tzz = fourier_vertical_down(gravity_tz)

    # Compute a source-valid stencil for every derivative cell, including the
    # finite-difference halo. Fourier fields are nearest-filled outside support,
    # but windows at real-data holes / the footprint edge are never accepted.
    raw_common = footprint & tmi_mask & tmi_vg_mask & grav_mask
    window_valid = raw_common
    mag_cloud = solve_euler_cloud(magnetic_tx, magnetic_ty, magnetic_tz, window_valid)
    grav_cloud = solve_euler_cloud(gravity_tx, gravity_ty, gravity_tzz, window_valid)
    pairs = pair_clouds(mag_cloud, grav_cloud)
    kde = kde_from_pairs(pairs, footprint.shape)
    pred, emission = emit_top_budget(kde, footprint, labels)
    pred[labels == 1] = 0.0
    pred[~footprint] = 0.0
    summary: dict[str, Any] = {
        "hypothesis_id": "H1",
        "transform": list(transform)[:6],
        "shape": list(footprint.shape),
        "feature_bands": {name: {"band": band_indices[name], "band_name": name} for name in ("tmi", "tmi_vg", "iso_grav_anom")},
        "magnetic_vertical_derivative_check": derivative_check,
        "magnetic_vertical_derivative_source": magnetic_vertical_source,
        "gravity_vertical_derivative_source": "padded reflect +|k| Fourier derivative of iso_grav_anom; supplied iso_grav_anom_vg not used",
        "common_real_input_pixels": int(common_valid.sum()),
        "exact_known_label_pixels_suppressed": int(np.count_nonzero((labels == 1) & footprint)),
        "magnetic_cloud": mag_cloud.to_summary(),
        "gravity_cloud": grav_cloud.to_summary(),
        "paired_cloud": pairs.to_summary(),
        "emission": emission,
    }
    return pred, summary
