"""H4: depth-aware dual-potential-field Euler solution-cloud KDE.

The Euler equation follows Reid et al. (1990), equation 5. For a measured field T at z=0,

    (x - x0) Tx + (y - y0) Ty - z0 Tz = N (B - T)

where z0 is positive down. For N=0 the regional term disappears and the active solve has three
unknowns. For N != 0, N*B is fitted as the fourth linear coefficient. H4 fixes N=0 for reduced-to-
pole magnetics and N=-1 for the raw gravity anomaly, with the latter explicitly treated as a
single-window approximation to the finite-step gravity case (not the generalized finite-step
inversion discussed by Reid & Thurston, 2014).

Unlike a gradient-threshold detector, this module returns source coordinates, depth, fit residual,
condition number, and the matched 3-D point cloud used to form the output KDE.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

from .research_euler import (
    RESOLUTION_M,
    _gradient_components,
    fill_nearest,
    fourier_vertical_down,
)


WINDOW_PX = 10
STRIDE_PX = 4
MAX_CONDITION = 1e5
MAX_RELATIVE_RESIDUAL = 0.20
MIN_DEPTH_M = 100.0
MAX_DEPTH_M = 2500.0
MAX_OFFSET_PX = 5.0
VOXEL_XY_M = 200.0
VOXEL_DEPTH_M = 200.0
LOCAL_XY_M = 300.0
LOCAL_DEPTH_M = 300.0
MIN_LOCAL_NEIGHBOURS = 3
DEPTH_MAD_SCALE_M = 150.0
SHALLOW_DECAY_M = 900.0
FIT_OFFSET_SCALE_PX = 2.0
PAIR_XY_M = 300.0
PAIR_DEPTH_M = 300.0
KDE_SIGMA_PX = 1.5
KDE_TRUNCATE = 3.0


@dataclass
class EulerCloud:
    """Accepted source solutions, all coordinates in grid pixels except depth."""

    row: np.ndarray
    col: np.ndarray
    depth_m: np.ndarray
    relative_residual: np.ndarray
    condition: np.ndarray
    offset_px: np.ndarray
    si: float
    field: str

    def __len__(self) -> int:
        return int(self.depth_m.size)

    def summary(self) -> dict[str, Any]:
        if not len(self):
            return {
                "field": self.field,
                "si": self.si,
                "accepted_solutions": 0,
            }
        return {
            "field": self.field,
            "si": self.si,
            "accepted_solutions": len(self),
            "depth_m_p05_p50_p95": np.percentile(self.depth_m, [5, 50, 95]).tolist(),
            "relative_residual_p50_p95": np.percentile(self.relative_residual, [50, 95]).tolist(),
            "condition_p50_p95": np.percentile(self.condition, [50, 95]).tolist(),
            "source_offset_px_p50_p95": np.percentile(self.offset_px, [50, 95]).tolist(),
        }


@dataclass
class VoxelCloud:
    """One representative solution per 3-D source voxel, plus its cluster weight."""

    row: np.ndarray
    col: np.ndarray
    x_m: np.ndarray
    y_m: np.ndarray
    depth_m: np.ndarray
    base_weight: np.ndarray
    cluster_weight: np.ndarray
    weight: np.ndarray
    neighbour_count: np.ndarray
    depth_mad_m: np.ndarray

    def __len__(self) -> int:
        return int(self.depth_m.size)


@dataclass
class PairCloud:
    row: np.ndarray
    col: np.ndarray
    depth_m: np.ndarray
    weight: np.ndarray
    xy_distance_m: np.ndarray
    depth_delta_m: np.ndarray

    def summary(self) -> dict[str, Any]:
        if not self.row.size:
            return {"mutual_nearest_pairs": 0, "weight_sum": 0.0}
        return {
            "mutual_nearest_pairs": int(self.row.size),
            "weight_sum": float(self.weight.sum(dtype=np.float64)),
            "xy_distance_m_p50_p95": np.percentile(self.xy_distance_m, [50, 95]).tolist(),
            "depth_delta_m_p50_p95": np.percentile(self.depth_delta_m, [50, 95]).tolist(),
            "pair_depth_m_p05_p50_p95": np.percentile(self.depth_m, [5, 50, 95]).tolist(),
        }


def potential_field_derivatives(field: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return filled T, its valid-cell mask, horizontal derivatives and downward-positive Tz.

    Invalid cells are nearest-filled solely for the derivative transform. Source inversion later
    rejects every window touching an invalid input cell (or its one-cell derivative halo).
    """
    filled, valid = fill_nearest(np.asarray(field, dtype=np.float32))
    tx, ty = _gradient_components(filled, resolution_m=RESOLUTION_M)
    tz = fourier_vertical_down(filled, pad_cells=128, resolution_m=RESOLUTION_M)
    return filled, valid, tx, ty, tz


def _empty_cloud(si: float, field: str) -> EulerCloud:
    empty = np.empty(0, dtype=np.float32)
    return EulerCloud(empty, empty.copy(), empty.copy(), empty.copy(), empty.copy(), empty.copy(), si, field)


def solve_euler_cloud(
    field_values: np.ndarray,
    tx: np.ndarray,
    ty: np.ndarray,
    tz: np.ndarray,
    valid: np.ndarray,
    *,
    transform,
    si: float,
    field: str,
    window_px: int = WINDOW_PX,
    stride_px: int = STRIDE_PX,
    chunk_window_rows: int = 8,
    cell_m: float = RESOLUTION_M,
) -> tuple[EulerCloud, dict[str, int]]:
    """Solve a vectorized, local-coordinate Euler system in sliding windows.

    The horizontal coordinates are centred on each window before the normal equations are formed,
    avoiding the ill-conditioning caused by using million-metre UTM coordinates as regression
    columns. SI=0 uses exactly the three identifiable location unknowns; non-zero SI includes the
    SI*B background coefficient. Window validity is checked at full resolution, not just on the
    subsampled Euler lattice.
    """
    T, tx, ty, tz = (np.asarray(a, dtype=np.float32) for a in (field_values, tx, ty, tz))
    valid = np.asarray(valid, dtype=bool)
    if T.ndim != 2 or any(a.shape != T.shape for a in (tx, ty, tz, valid)):
        raise ValueError("field, derivatives, and valid mask must have identical 2-D shapes")
    if si not in (0.0, -1.0):
        raise ValueError("H4 locks magnetic SI=0 and gravity SI=-1")
    if window_px < 5 or stride_px <= 0 or chunk_window_rows <= 0:
        raise ValueError("window must be at least 5 px; stride and chunk rows must be positive")

    H, W = T.shape
    if H < window_px or W < window_px:
        return _empty_cloud(si, field), {"windows_total": 0, "windows_input_valid": 0, "windows_conditioned": 0, "solutions_accepted": 0}

    # A central-difference halo is required for every input pixel used by a window.
    derivative_valid = ndimage.minimum_filter(valid.astype(np.uint8), size=3, mode="constant", cval=0) > 0
    valid_windows = np.lib.stride_tricks.sliding_window_view(
        derivative_valid, (window_px, window_px)
    )[::stride_px, ::stride_px].all(axis=(-1, -2))
    txw = np.lib.stride_tricks.sliding_window_view(tx, (window_px, window_px))[::stride_px, ::stride_px]
    tyw = np.lib.stride_tricks.sliding_window_view(ty, (window_px, window_px))[::stride_px, ::stride_px]
    tzw = np.lib.stride_tricks.sliding_window_view(tz, (window_px, window_px))[::stride_px, ::stride_px]
    tw = np.lib.stride_tricks.sliding_window_view(T, (window_px, window_px))[::stride_px, ::stride_px]
    n_rows, n_cols = valid_windows.shape

    xoff = (np.arange(window_px, dtype=np.float32) + 0.5 - window_px / 2.0) * cell_m
    xgrid = np.broadcast_to(xoff[None, :], (window_px, window_px))
    ygrid = np.broadcast_to((-xoff)[:, None], (window_px, window_px))
    n_unknowns = 3 if si == 0.0 else 4

    row_parts: list[np.ndarray] = []
    col_parts: list[np.ndarray] = []
    depth_parts: list[np.ndarray] = []
    residual_parts: list[np.ndarray] = []
    condition_parts: list[np.ndarray] = []
    offset_parts: list[np.ndarray] = []
    tested = int(valid_windows.size)
    input_valid = int(valid_windows.sum())
    conditioned = 0
    solved_finite = 0
    depth_gate_count = 0
    fit_gate_count = 0
    offset_gate_count = 0
    all_physical_gates_count = 0

    for r0 in range(0, n_rows, chunk_window_rows):
        r1 = min(n_rows, r0 + chunk_window_rows)
        ax, ay, az, at = txw[r0:r1], tyw[r0:r1], tzw[r0:r1], tw[r0:r1]
        rhs_image = xgrid * ax + ygrid * ay
        if si != 0.0:
            rhs_image = rhs_image + np.float32(si) * at
        rhs = rhs_image.reshape(r1 - r0, n_cols, window_px * window_px).astype(np.float64)
        rhs_norm = np.einsum("...k,...k->...", rhs, rhs, optimize=True)

        if si == 0.0:
            design = np.stack((ax, ay, az), axis=-1).reshape(
                r1 - r0, n_cols, window_px * window_px, 3
            ).astype(np.float64)
        else:
            intercept = np.full(ax.shape, np.float32(si), dtype=np.float32)
            design = np.stack((ax, ay, az, intercept), axis=-1).reshape(
                r1 - r0, n_cols, window_px * window_px, 4
            ).astype(np.float64)
        normal = np.einsum("...ki,...kj->...ij", design, design, optimize=True)
        cross = np.einsum("...ki,...k->...i", design, rhs, optimize=True)
        eigenvalues = np.linalg.eigvalsh(normal)
        condition = np.full(eigenvalues.shape[:-1], np.inf, dtype=np.float64)
        eig_ok = eigenvalues[..., 0] > 0.0
        condition[eig_ok] = eigenvalues[..., -1][eig_ok] / eigenvalues[..., 0][eig_ok]
        finite = (
            np.isfinite(rhs_norm)
            & np.isfinite(condition)
            & (condition <= MAX_CONDITION)
            & (rhs_norm > 1e-24)
        )
        conditioned += int(np.count_nonzero(finite & valid_windows[r0:r1]))
        keep_window = valid_windows[r0:r1] & finite
        if not keep_window.any():
            continue

        normal_good = normal[keep_window]
        cross_good = cross[keep_window]
        try:
            solutions = np.linalg.solve(normal_good, cross_good[..., None])[..., 0]
        except np.linalg.LinAlgError:
            solved: list[np.ndarray] = []
            for matrix, vector in zip(normal_good, cross_good):
                try:
                    solved.append(np.linalg.solve(matrix, vector))
                except np.linalg.LinAlgError:
                    solved.append(np.full(n_unknowns, np.nan, dtype=np.float64))
            solutions = np.asarray(solved, dtype=np.float64)

        fit = np.einsum("...i,...i->...", cross_good, solutions, optimize=True)
        residual_sq = np.maximum(rhs_norm[keep_window] - fit, 0.0)
        residual = np.sqrt(residual_sq / (rhs_norm[keep_window] + 1e-30))
        loc = np.argwhere(keep_window)
        starts_row = (loc[:, 0] + r0) * stride_px
        starts_col = loc[:, 1] * stride_px
        dx0, dy0, depth = solutions[:, 0], solutions[:, 1], solutions[:, 2]
        offset_px = np.hypot(dx0, dy0) / cell_m
        finite_solution = np.isfinite(dx0) & np.isfinite(dy0) & np.isfinite(depth) & np.isfinite(residual)
        depth_gate = finite_solution & (depth >= MIN_DEPTH_M) & (depth <= MAX_DEPTH_M)
        fit_gate = finite_solution & (residual <= MAX_RELATIVE_RESIDUAL)
        offset_gate = finite_solution & (offset_px <= MAX_OFFSET_PX)
        accepted = (
            finite_solution
            & np.isfinite(condition[keep_window])
            & (depth >= MIN_DEPTH_M)
            & (depth <= MAX_DEPTH_M)
            & (residual <= MAX_RELATIVE_RESIDUAL)
            & (offset_px <= MAX_OFFSET_PX)
        )
        solved_finite += int(finite_solution.sum())
        depth_gate_count += int(depth_gate.sum())
        fit_gate_count += int(fit_gate.sum())
        offset_gate_count += int(offset_gate.sum())
        all_physical_gates_count += int(accepted.sum())
        if not accepted.any():
            continue
        centre_offset = (window_px - 1.0) / 2.0
        row_parts.append((starts_row[accepted] + centre_offset - dy0[accepted] / cell_m).astype(np.float32))
        col_parts.append((starts_col[accepted] + centre_offset + dx0[accepted] / cell_m).astype(np.float32))
        depth_parts.append(depth[accepted].astype(np.float32))
        residual_parts.append(residual[accepted].astype(np.float32))
        condition_parts.append(condition[keep_window][accepted].astype(np.float32))
        offset_parts.append(offset_px[accepted].astype(np.float32))

    if row_parts:
        cloud = EulerCloud(
            row=np.concatenate(row_parts),
            col=np.concatenate(col_parts),
            depth_m=np.concatenate(depth_parts),
            relative_residual=np.concatenate(residual_parts),
            condition=np.concatenate(condition_parts),
            offset_px=np.concatenate(offset_parts),
            si=si,
            field=field,
        )
    else:
        cloud = _empty_cloud(si, field)
    counts = {
        "windows_total": tested,
        "windows_input_valid": input_valid,
        "windows_conditioned": conditioned,
        "solutions_finite_before_physical_gates": solved_finite,
        "solutions_depth_gate": depth_gate_count,
        "solutions_relative_fit_gate": fit_gate_count,
        "solutions_offset_gate": offset_gate_count,
        "solutions_pass_all_physical_gates": all_physical_gates_count,
        "solutions_accepted": len(cloud),
    }
    return cloud, counts


def _filter_to_footprint(cloud: EulerCloud, footprint: np.ndarray) -> np.ndarray:
    h, w = footprint.shape
    r = np.floor(cloud.row.astype(np.float64) + 0.5).astype(np.int64)
    c = np.floor(cloud.col.astype(np.float64) + 0.5).astype(np.int64)
    inside = (r >= 0) & (r < h) & (c >= 0) & (c < w)
    keep = np.zeros(len(cloud), dtype=bool)
    keep[inside] = footprint[r[inside], c[inside]]
    return keep


def voxelize_and_weight(
    cloud: EulerCloud,
    footprint: np.ndarray,
    transform,
    *,
    voxel_xy_m: float = VOXEL_XY_M,
    voxel_depth_m: float = VOXEL_DEPTH_M,
    local_xy_m: float = LOCAL_XY_M,
    local_depth_m: float = LOCAL_DEPTH_M,
    min_neighbours: int = MIN_LOCAL_NEIGHBOURS,
    depth_mad_scale_m: float = DEPTH_MAD_SCALE_M,
    shallow_decay_m: float = SHALLOW_DECAY_M,
    fit_offset_scale_px: float = FIT_OFFSET_SCALE_PX,
    source_context_factor: np.ndarray | None = None,
) -> tuple[VoxelCloud, dict[str, Any]]:
    """Deduplicate overlap-votes and weight compact local 3-D clusters."""
    empty = np.empty(0, dtype=np.float32)
    empty_cloud = VoxelCloud(*(empty.copy() for _ in range(10)))
    if not len(cloud):
        return empty_cloud, {"raw_solutions": 0, "in_footprint_solutions": 0, "unique_source_voxels": 0, "cluster_supported_voxels": 0}

    footprint_keep = _filter_to_footprint(cloud, footprint)
    idx = np.flatnonzero(footprint_keep)
    if idx.size == 0:
        return empty_cloud, {"raw_solutions": len(cloud), "in_footprint_solutions": 0, "unique_source_voxels": 0, "cluster_supported_voxels": 0}

    if source_context_factor is None:
        context = np.ones(len(cloud), dtype=np.float64)
    else:
        context = np.asarray(source_context_factor, dtype=np.float64)
        if context.shape != (len(cloud),):
            raise ValueError("source_context_factor must have one value per Euler solution")
        if not np.isfinite(context).all() or (context < 0.0).any() or (context > 1.0).any():
            raise ValueError("source_context_factor values must be finite and in [0,1]")
    row = cloud.row[idx].astype(np.float64)
    col = cloud.col[idx].astype(np.float64)
    depth = cloud.depth_m[idx].astype(np.float64)
    x_m = transform.c + (col + 0.5) * transform.a
    y_m = transform.f + (row + 0.5) * transform.e
    residual = cloud.relative_residual[idx].astype(np.float64)
    offset = cloud.offset_px[idx].astype(np.float64)
    context = context[idx]
    base = (
        np.exp(-depth / shallow_decay_m)
        * np.exp(-0.5 * (residual / MAX_RELATIVE_RESIDUAL) ** 2)
        * np.exp(-0.5 * (offset / fit_offset_scale_px) ** 2)
        * context
    )

    bx = np.floor(x_m / voxel_xy_m).astype(np.int64)
    by = np.floor(y_m / voxel_xy_m).astype(np.int64)
    bz = np.floor(depth / voxel_depth_m).astype(np.int64)
    original_order = np.arange(idx.size, dtype=np.int64)
    order = np.lexsort((original_order, -base, bz, by, bx))
    sx, sy, sz = bx[order], by[order], bz[order]
    first = np.r_[True, (sx[1:] != sx[:-1]) | (sy[1:] != sy[:-1]) | (sz[1:] != sz[:-1])]
    selected = order[first]

    row = row[selected].astype(np.float32)
    col = col[selected].astype(np.float32)
    depth = depth[selected].astype(np.float32)
    x_m = x_m[selected].astype(np.float64)
    y_m = y_m[selected].astype(np.float64)
    base = base[selected].astype(np.float32)

    if depth.size < min_neighbours:
        supported = np.zeros(depth.size, dtype=bool)
        neighbour_count = np.zeros(depth.size, dtype=np.int16)
        depth_mad = np.full(depth.size, np.inf, dtype=np.float32)
        cluster_weight = np.zeros(depth.size, dtype=np.float32)
    else:
        points = np.column_stack((x_m, y_m, depth.astype(np.float64)))
        tree = cKDTree(points)
        radius = float(np.hypot(local_xy_m, local_depth_m))
        neighbours = tree.query_ball_point(points, r=radius, workers=-1)
        neighbour_count = np.zeros(depth.size, dtype=np.int16)
        depth_mad = np.full(depth.size, np.inf, dtype=np.float32)
        for i, candidates in enumerate(neighbours):
            cand = np.asarray(candidates, dtype=np.int64)
            dx = points[cand, 0] - points[i, 0]
            dy = points[cand, 1] - points[i, 1]
            dz = points[cand, 2] - points[i, 2]
            local = (np.hypot(dx, dy) <= local_xy_m) & (np.abs(dz) <= local_depth_m)
            cand = cand[local]
            neighbour_count[i] = min(int(cand.size), np.iinfo(np.int16).max)
            if cand.size:
                d = depth[cand].astype(np.float64)
                med = float(np.median(d))
                depth_mad[i] = float(np.median(np.abs(d - med)))
        supported = neighbour_count >= min_neighbours
        density = np.minimum(neighbour_count.astype(np.float32) / float(min_neighbours), 1.0)
        cluster_weight = density * np.exp(-np.nan_to_num(depth_mad, nan=np.inf, posinf=np.inf) / depth_mad_scale_m)
        cluster_weight[~supported] = 0.0

    total_weight = base * cluster_weight
    voxel_cloud = VoxelCloud(
        row=row,
        col=col,
        x_m=x_m,
        y_m=y_m,
        depth_m=depth,
        base_weight=base,
        cluster_weight=cluster_weight.astype(np.float32),
        weight=total_weight.astype(np.float32),
        neighbour_count=neighbour_count,
        depth_mad_m=depth_mad.astype(np.float32),
    )
    stats = {
        "raw_solutions": int(len(cloud)),
        "in_footprint_solutions": int(idx.size),
        "unique_source_voxels": int(len(voxel_cloud)),
        "cluster_supported_voxels": int(np.count_nonzero(total_weight > 0)),
        "neighbour_count_p50_p95": np.percentile(neighbour_count, [50, 95]).tolist() if neighbour_count.size else [],
        "depth_mad_m_p50_p95_supported": np.percentile(depth_mad[supported], [50, 95]).tolist() if supported.any() else [],
        "source_context_factor_p05_p50_p95": np.percentile(context, [5, 50, 95]).tolist() if context.size else [],
    }
    return voxel_cloud, stats


def mutual_nearest_pairs(
    magnetic: VoxelCloud,
    gravity: VoxelCloud,
    *,
    max_xy_m: float = PAIR_XY_M,
    max_depth_m: float = PAIR_DEPTH_M,
) -> PairCloud:
    """Return one-to-one mutual-nearest 3-D magnetic/gravity source pairs."""
    if len(magnetic) == 0 or len(gravity) == 0:
        empty = np.empty(0, dtype=np.float32)
        return PairCloud(*(empty.copy() for _ in range(6)))
    mi = np.flatnonzero(magnetic.weight > 0)
    gi = np.flatnonzero(gravity.weight > 0)
    if mi.size == 0 or gi.size == 0:
        empty = np.empty(0, dtype=np.float32)
        return PairCloud(*(empty.copy() for _ in range(6)))

    m_xyz = np.column_stack((magnetic.x_m[mi] / max_xy_m, magnetic.y_m[mi] / max_xy_m, magnetic.depth_m[mi] / max_depth_m))
    g_xyz = np.column_stack((gravity.x_m[gi] / max_xy_m, gravity.y_m[gi] / max_xy_m, gravity.depth_m[gi] / max_depth_m))
    mag_tree = cKDTree(m_xyz)
    grav_tree = cKDTree(g_xyz)
    _, nearest_g = grav_tree.query(m_xyz, k=1, workers=-1)
    _, nearest_m = mag_tree.query(g_xyz, k=1, workers=-1)
    m_local = np.arange(mi.size, dtype=np.int64)
    reciprocal = nearest_m[nearest_g] == m_local
    gm_local = nearest_g
    m_sel = m_local[reciprocal]
    g_sel = gm_local[reciprocal]
    if m_sel.size == 0:
        empty = np.empty(0, dtype=np.float32)
        return PairCloud(*(empty.copy() for _ in range(6)))

    m_idx, g_idx = mi[m_sel], gi[g_sel]
    dx = magnetic.x_m[m_idx] - gravity.x_m[g_idx]
    dy = magnetic.y_m[m_idx] - gravity.y_m[g_idx]
    dz = magnetic.depth_m[m_idx] - gravity.depth_m[g_idx]
    dxy = np.hypot(dx, dy)
    keep = (dxy <= max_xy_m) & (np.abs(dz) <= max_depth_m)
    m_idx, g_idx = m_idx[keep], g_idx[keep]
    dxy, dz = dxy[keep], dz[keep]
    if m_idx.size == 0:
        empty = np.empty(0, dtype=np.float32)
        return PairCloud(*(empty.copy() for _ in range(6)))

    penalty = np.exp(-0.5 * ((dxy / max_xy_m) ** 2 + (dz / max_depth_m) ** 2))
    pair_weight = np.sqrt(magnetic.weight[m_idx].astype(np.float64) * gravity.weight[g_idx].astype(np.float64)) * penalty
    row = 0.5 * (magnetic.row[m_idx].astype(np.float64) + gravity.row[g_idx].astype(np.float64))
    col = 0.5 * (magnetic.col[m_idx].astype(np.float64) + gravity.col[g_idx].astype(np.float64))
    depth = 0.5 * (magnetic.depth_m[m_idx].astype(np.float64) + gravity.depth_m[g_idx].astype(np.float64))
    return PairCloud(
        row=row.astype(np.float32),
        col=col.astype(np.float32),
        depth_m=depth.astype(np.float32),
        weight=pair_weight.astype(np.float32),
        xy_distance_m=dxy.astype(np.float32),
        depth_delta_m=np.abs(dz).astype(np.float32),
    )


def _bilinear_deposit(
    shape: tuple[int, int],
    footprint: np.ndarray,
    row: np.ndarray,
    col: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    acc = np.zeros(shape, dtype=np.float64)
    r0 = np.floor(row).astype(np.int64)
    c0 = np.floor(col).astype(np.int64)
    dr = row - r0
    dc = col - c0
    for rr, cc, ww in (
        (r0, c0, (1.0 - dr) * (1.0 - dc)),
        (r0 + 1, c0, dr * (1.0 - dc)),
        (r0, c0 + 1, (1.0 - dr) * dc),
        (r0 + 1, c0 + 1, dr * dc),
    ):
        inside = (rr >= 0) & (rr < shape[0]) & (cc >= 0) & (cc < shape[1])
        loc = np.flatnonzero(inside)
        if loc.size:
            loc = loc[footprint[rr[loc], cc[loc]]]
            np.add.at(acc, (rr[loc], cc[loc]), (weights * ww)[loc])
    return acc


def project_weighted_cloud_to_kde(
    row: np.ndarray,
    col: np.ndarray,
    depth_m: np.ndarray,
    weight: np.ndarray,
    footprint: np.ndarray,
    known_labels: np.ndarray,
    *,
    sigma_px: float = KDE_SIGMA_PX,
    truncate: float = KDE_TRUNCATE,
    stop_message: str = "H4 stop: no positive 3-D source support",
) -> tuple[np.ndarray, dict[str, Any]]:
    """Project any depth-weighted 3-D source cloud through a compact continuous 2-D KDE."""
    footprint = np.asarray(footprint, dtype=bool)
    known = np.asarray(known_labels) == 1
    row = np.asarray(row, dtype=np.float32)
    col = np.asarray(col, dtype=np.float32)
    depth_m = np.asarray(depth_m, dtype=np.float32)
    weight = np.asarray(weight, dtype=np.float32)
    if footprint.shape != known.shape:
        raise ValueError("footprint and known-label mask shapes differ")
    if any(x.shape != row.shape for x in (col, depth_m, weight)):
        raise ValueError("source row, column, depth, and weight arrays must have identical shapes")
    if sigma_px <= 0 or truncate <= 0:
        raise ValueError("KDE sigma and truncate must be positive")
    if not row.size or not np.any(weight > 0):
        raise RuntimeError(stop_message)
    if not (np.isfinite(row).all() and np.isfinite(col).all() and np.isfinite(depth_m).all() and np.isfinite(weight).all()):
        raise ValueError("positive source cloud coordinates and weights must be finite")
    acc = _bilinear_deposit(footprint.shape, footprint, row, col, weight)
    kde = ndimage.gaussian_filter(acc, sigma=sigma_px, mode="constant", truncate=truncate)
    kde[~footprint] = 0.0
    kde[known & footprint] = 0.0
    maximum = float(kde[footprint].max()) if footprint.any() else 0.0
    if not np.isfinite(maximum) or maximum <= 0.0:
        raise RuntimeError("source KDE has no positive in-footprint support after exact label suppression")
    prediction = (kde / maximum).astype(np.float32)
    prediction[~footprint] = np.nan
    prediction[known & footprint] = 0.0
    return prediction, {
        "source_points": int(row.size),
        "source_weight_sum": float(weight.sum(dtype=np.float64)),
        "source_depth_m_p05_p50_p95": np.percentile(depth_m, [5, 50, 95]).tolist(),
        "kde_sigma_px": float(sigma_px),
        "kde_truncate_sigma": float(truncate),
        "normalization": "divide by maximum in-footprint KDE after exact-label suppression",
        "positive_in_footprint": int(np.count_nonzero(prediction[footprint] > 0)),
        "max_in_footprint": float(np.max(prediction[footprint])),
        "known_label_pixels_forced_zero": int(np.count_nonzero(prediction[known & footprint] == 0)),
        "outside_footprint_nan_required_on_write": int((~footprint).sum()),
    }


def project_pair_cloud_to_kde(
    pairs: PairCloud,
    footprint: np.ndarray,
    known_labels: np.ndarray,
    *,
    sigma_px: float = KDE_SIGMA_PX,
    truncate: float = KDE_TRUNCATE,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Bilinearly deposit 3-D matched source pairs and project through a compact 2-D KDE."""
    if not pairs.row.size or not np.any(pairs.weight > 0):
        raise RuntimeError("H4 stop: locked Euler/pair gates produced no positive 3-D pair support")
    prediction, info = project_weighted_cloud_to_kde(
        pairs.row, pairs.col, pairs.depth_m, pairs.weight, footprint, known_labels,
        sigma_px=sigma_px, truncate=truncate,
    )
    return prediction, {
        **info,
        "pair_count": int(pairs.row.size),
        "pair_weight_sum": float(pairs.weight.sum(dtype=np.float64)),
    }


def build_h4_prediction(
    magnetic: EulerCloud,
    gravity: EulerCloud,
    footprint: np.ndarray,
    known_labels: np.ndarray,
    transform,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Apply locked H4 voxel weighting, mutual 3-D matching, and continuous KDE projection."""
    mag_voxels, mag_stats = voxelize_and_weight(magnetic, footprint, transform)
    grav_voxels, grav_stats = voxelize_and_weight(gravity, footprint, transform)
    pairs = mutual_nearest_pairs(mag_voxels, grav_voxels)
    prediction, kde_stats = project_pair_cloud_to_kde(pairs, footprint, known_labels)
    summary = {
        "hypothesis": "H4",
        "construction_only": True,
        "not_scored_during_construction": True,
        "magnetic": {**magnetic.summary(), "voxelization": mag_stats},
        "gravity": {**gravity.summary(), "voxelization": grav_stats},
        "mutual_nearest_pairing": pairs.summary(),
        "kde": kde_stats,
        "locked_parameters": {
            "magnetic_layer": "rtp",
            "magnetic_si": 0,
            "gravity_layer": "iso_grav_anom",
            "gravity_si": -1,
            "window_px": WINDOW_PX,
            "stride_px": STRIDE_PX,
            "max_condition": MAX_CONDITION,
            "max_relative_residual": MAX_RELATIVE_RESIDUAL,
            "depth_range_m": [MIN_DEPTH_M, MAX_DEPTH_M],
            "max_source_offset_px": MAX_OFFSET_PX,
            "voxel_xy_m": VOXEL_XY_M,
            "voxel_depth_m": VOXEL_DEPTH_M,
            "local_neighbour_xy_m": LOCAL_XY_M,
            "local_neighbour_depth_m": LOCAL_DEPTH_M,
            "minimum_neighbours_including_self": MIN_LOCAL_NEIGHBOURS,
            "depth_mad_scale_m": DEPTH_MAD_SCALE_M,
            "shallow_decay_m": SHALLOW_DECAY_M,
            "pair_xy_tolerance_m": PAIR_XY_M,
            "pair_depth_tolerance_m": PAIR_DEPTH_M,
            "kde_sigma_px": KDE_SIGMA_PX,
            "kde_truncate_sigma": KDE_TRUNCATE,
            "known_label_rule": "zero exact known-label pixels only; no distance buffer",
        },
    }
    return prediction, summary


def summarize_inputs_and_derivatives(field: np.ndarray, valid: np.ndarray, name: str) -> dict[str, Any]:
    values = np.asarray(field)[valid]
    return {
        "field": name,
        "valid_cells": int(valid.sum()),
        "field_min_max": [float(values.min()), float(values.max())] if values.size else [None, None],
        "field_std": float(values.std(dtype=np.float64)) if values.size else None,
    }
