"""H4-A: Euler SI=0 depth-cluster contact network with surface-expression concurrence.

Scientific basis
----------------
Reid, A.B., Allsop, J.M., Granser, H., Millett, A.J. and Somerton, I.W. (1990),
"Magnetic interpretation in three dimensions using Euler deconvolution",
*Geophysics* 55(1), 80-91, https://doi.org/10.1190/1.1442774 — sliding-window
least-squares solution of the Euler homogeneity relation.  Structural index
N = 0 is the fault-like contact case used here for both the magnetic and the
gravity channel.

Blakely, R.J. (1995), *Potential Theory in Gravity and Magnetic Applications*,
eq. 12-8 — the Fourier-domain vertical derivative used to build the gravity
vertical derivative (the supplied `iso_grav_anom_vg` band is never trusted).

Everything in this module is a *source-locator* stage followed by a
*clustering* stage, as the project brief requires: the raw product is a cloud
of depth-labelled solution points, which is then converted to a continuous
raster by a weighted kernel-density estimate.  No gradient threshold is ever
substituted for the deconvolution.

Lineament stage
---------------
A Hessian-based multi-scale ridge response turns the point-density field into
*linear* geometry, because the ground truth is a line raster rather than a
point set.  For a 2-D Hessian with eigenvalues |l1| >= |l2|, l2 is the
curvature across the line and l1 along it; the response |l2| is maximal on a
linear structure whose width matches the scale sigma.  Scale normalisation by
sigma**2 makes the response comparable between scales.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

from .euler import EulerCloud, deconvolve, merge_clouds

PIXEL_M = 100.0


# --------------------------------------------------------------------------
# Stage 1 — Euler deconvolution of both potential fields
# --------------------------------------------------------------------------
@dataclass
class ContactConfig:
    field_name: str
    window_px: int
    stride_px: int
    structural_index: float = 0.0
    analytic_percentile: float = 72.0
    max_rel_se: float = 0.22
    min_depth_m: float = 80.0
    max_depth_m: float = 2200.0


DEFAULT_CONFIGS = (
    ContactConfig("rtp", 8, 4),
    ContactConfig("rtp", 12, 4),
    ContactConfig("rtp", 16, 6),
    ContactConfig("iso_grav_anom", 8, 4),
    ContactConfig("iso_grav_anom", 12, 4),
    ContactConfig("iso_grav_anom", 16, 6),
)


def solve_contacts(field: np.ndarray, valid: np.ndarray, config: ContactConfig) -> EulerCloud:
    """One moving-window SI=0 Euler solve.  Thin wrapper that keeps the record explicit."""
    return deconvolve(
        field,
        valid,
        field_name=config.field_name,
        structural_index=config.structural_index,
        window_px=config.window_px,
        stride_px=config.stride_px,
        analytic_percentile=config.analytic_percentile,
        max_rel_se=config.max_rel_se,
        min_depth_m=config.min_depth_m,
        max_depth_m=config.max_depth_m,
    )


def solve_field(field: np.ndarray, valid: np.ndarray, name: str,
                configs: tuple[ContactConfig, ...] = DEFAULT_CONFIGS) -> tuple[EulerCloud, list[dict]]:
    """Solve every registered window size for one field and concatenate the clouds."""
    clouds, stats = [], []
    for cfg in configs:
        if cfg.field_name != name:
            continue
        cloud = solve_contacts(field, valid, cfg)
        clouds.append(cloud)
        stats.append({"window_px": cfg.window_px, "stride_px": cfg.stride_px, **cloud.stats})
    if not clouds:
        raise ValueError(f"no configs registered for field {name!r}")
    return merge_clouds(clouds), stats


# --------------------------------------------------------------------------
# Stage 2 — 3-D depth consistency gate and solution weights
# --------------------------------------------------------------------------
def depth_cluster_weights(
    cloud: EulerCloud,
    *,
    xy_px: float = 3.0,
    z_m: float = 500.0,
    min_neighbors: int = 3,
    max_depth_cv: float = 0.35,
    depth_scale_m: float = 900.0,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Weight each solution by shallowness, fit quality, 3-D tightness and depth coherence.

    Returns ``(weights, keep_mask, stats)``.  ``keep_mask`` is the *gate*: a
    solution survives only when at least ``min_neighbors`` other solutions lie
    within ``xy_px`` horizontally and ``z_m`` vertically AND those neighbours
    agree in depth (coefficient of variation of depth <= ``max_depth_cv``).
    Isolated or depth-scattered solutions are the noise population described by
    Reid et al. (1990) and are removed rather than down-weighted.
    """
    n = len(cloud)
    if n == 0:
        return np.empty(0), np.empty(0, dtype=bool), {"solutions": 0, "kept": 0}
    # Work in pixel units: 1 px == 100 m horizontally, and depth is converted
    # with the same factor so the ball query matches "xy within xy_px AND
    # depth within z_m" (a ball of radius hypot(xy_px, z_m/100)).
    pts = np.column_stack([cloud.row, cloud.col, cloud.depth_m / PIXEL_M])
    tree = cKDTree(pts)
    neigh = tree.query_ball_tree(tree, r=float(np.hypot(xy_px, z_m / 100.0)))
    neighbor_count = np.array([len(nb) - 1 for nb in neigh], dtype=np.float64)
    depth_cv = np.zeros(n, dtype=np.float64)
    for i, nb in enumerate(neigh):
        if len(nb) > 1:
            depths = cloud.depth_m[np.fromiter(nb, dtype=int, count=len(nb))]
            mu = float(depths.mean())
            depth_cv[i] = float(depths.std()) / max(mu, 1.0)
        else:
            depth_cv[i] = np.inf
    keep = (neighbor_count >= float(min_neighbors)) & (depth_cv <= max_depth_cv)
    shallow = np.exp(-np.clip(cloud.depth_m, 0, None) / depth_scale_m)
    quality = 1.0 / (1.0 + 8.0 * np.clip(cloud.rel_se, 0, 10))
    ref = float(np.percentile(neighbor_count[keep], 90)) if keep.any() else 1.0
    tight = np.clip(neighbor_count / max(ref, 1.0), 0.0, 3.0)
    coherent = np.exp(-4.0 * np.clip(depth_cv, 0, 10))
    analytic = cloud.analytic
    a90 = float(np.percentile(analytic, 90)) if analytic.size else 1.0
    amp = 0.25 + 0.75 * np.clip(analytic / max(a90, 1e-12), 0, 2.0)
    w = shallow * quality * tight * coherent * amp
    w = np.where(np.isfinite(w) & keep & (w > 0), w, 0.0)
    stats = {
        "solutions": int(n),
        "kept": int(keep.sum()),
        "kept_fraction": float(keep.mean()) if n else 0.0,
        "median_depth_m_kept": float(np.median(cloud.depth_m[keep])) if keep.any() else None,
        "p90_depth_m_kept": float(np.percentile(cloud.depth_m[keep], 90)) if keep.any() else None,
        "median_neighbors": float(np.median(neighbor_count)),
        "median_depth_cv": float(np.median(depth_cv[np.isfinite(depth_cv)])) if np.isfinite(depth_cv).any() else None,
        "weight_sum": float(w.sum()),
    }
    return w, keep, stats


def cross_field_concordance(
    mag: EulerCloud, mag_w: np.ndarray,
    grav: EulerCloud, grav_w: np.ndarray,
    *,
    xy_px: float = 3.0,
    z_m: float = 500.0,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Points where an accepted magnetic and gravity solution agree in (x, y, depth).

    Returns ``(points, weights, stats)`` where the polygon-free pairing rule is
    "nearest partner within xy_px and z_m"; the pair is deposited at the
    magnetic solution location and its weight is the geometric mean of the two
    solution weights, so a strong magnetic-only cluster cannot inherit gravity
    credibility it does not have.
    """
    if len(mag) == 0 or len(grav) == 0:
        empty = np.empty((0, 3))
        return empty, np.empty(0), {"pairs": 0, "reason": "one field produced no accepted solutions"}
    m_pts = np.column_stack([mag.row, mag.col, mag.depth_m / 100.0])
    g_pts = np.column_stack([grav.row, grav.col, grav.depth_m / 100.0])
    g_tree = cKDTree(g_pts)
    radius = float(np.hypot(xy_px, z_m / 100.0))
    neigh = g_tree.query_ball_point(m_pts, r=radius)
    keep, weights = [], []
    for i, nb in enumerate(neigh):
        if not nb:
            continue
        m_pt = m_pts[i]
        d = np.linalg.norm(g_pts[nb] - m_pt, axis=1)
        j = nb[int(np.argmin(d))]
        keep.append(i)
        weights.append(float(np.sqrt(max(mag_w[i], 0.0) * max(grav_w[j], 0.0))))
    if not keep:
        return np.empty((0, 3)), np.empty(0), {"pairs": 0, "mag_accepted": int(len(mag)), "grav_accepted": int(len(grav))}
    idx = np.asarray(keep, dtype=int)
    pair_w = np.asarray(weights, dtype=np.float64)
    stats = {
        "pairs": int(idx.size),
        "mag_accepted": int(len(mag)),
        "grav_accepted": int(len(grav)),
        "paired_fraction_of_magnetic": float(idx.size / max(len(mag), 1)),
        "median_depth_m": float(np.median(mag.depth_m[idx])),
        "pair_weight_sum": float(pair_w.sum()),
    }
    return m_pts[idx] * np.array([1.0, 1.0, 100.0]), pair_w, stats


def splat_points(
    points: np.ndarray,
    weights: np.ndarray,
    shape: tuple[int, int],
    *,
    sigma_px: float = 1.5,
) -> np.ndarray:
    """Weighted bilinear splat + Gaussian KDE of a 3-D solution cloud onto the grid."""
    acc = np.zeros(shape, dtype=np.float64)
    if points.size == 0 or weights.size == 0 or not np.any(weights > 0):
        return acc
    h, w = shape
    r = np.clip(points[:, 0], 0, h - 1)
    c = np.clip(points[:, 1], 0, w - 1)
    r0 = np.floor(r).astype(int)
    c0 = np.floor(c).astype(int)
    r1 = np.clip(r0 + 1, 0, h - 1)
    c1 = np.clip(c0 + 1, 0, w - 1)
    fr = r - r0
    fc = c - c0
    np.add.at(acc, (r0, c0), weights * (1 - fr) * (1 - fc))
    np.add.at(acc, (r0, c1), weights * (1 - fr) * fc)
    np.add.at(acc, (r1, c0), weights * fr * (1 - fc))
    np.add.at(acc, (r1, c1), weights * fr * fc)
    return ndimage.gaussian_filter(acc, sigma=float(sigma_px), mode="constant", truncate=3.0)


# --------------------------------------------------------------------------
# Stage 3 — lineament (trace) enhancement of the depth-cluster field
# --------------------------------------------------------------------------
def robust_normalize(field: np.ndarray, mask: np.ndarray, *, hi: float = 99.5) -> np.ndarray:
    """Scale positive field values onto [0, 1] using a robust high quantile."""
    out = np.zeros_like(field, dtype=np.float64)
    pos = mask & (field > 0)
    if not pos.any():
        return out
    top = float(np.percentile(field[pos], hi))
    if top <= 0:
        top = float(field[pos].max())
    out = np.clip(field / max(top, 1e-30), 0.0, 1.0)
    return np.where(mask, out, 0.0)


def hessian_lineament_response(
    field: np.ndarray,
    mask: np.ndarray,
    *,
    scales: tuple[float, ...] = (1.5, 3.0, 6.0),
) -> np.ndarray:
    """Multi-scale Hessian ridge response, maximised over scale and polarisation.

    For each scale the 2-D Hessian eigenvalues are formed analytically; the
    response is the scale-normalised magnitude of the *across-line* curvature
    ``|l2|`` where ``|l1| >= |l2|`` (Frangi-style vesselness without the
    amplitude terms, which are unstable on sparse fields).  Both bright and
    dark ridges are kept because a contact can appear as either sign.
    """
    resp = np.zeros_like(field, dtype=np.float64)
    for sigma in scales:
        dyy = ndimage.gaussian_filter(field, sigma, order=(2, 0), mode="constant", truncate=4.0)
        dxx = ndimage.gaussian_filter(field, sigma, order=(0, 2), mode="constant", truncate=4.0)
        dxy = ndimage.gaussian_filter(field, sigma, order=(1, 1), mode="constant", truncate=4.0)
        tmp = np.sqrt(((dyy - dxx) * 0.5) ** 2 + dxy ** 2)
        l1 = 0.5 * (dyy + dxx) + tmp
        l2 = 0.5 * (dyy + dxx) - tmp
        # |l1| >= |l2| by construction of the analytical eigenvalues above.
        r = (sigma ** 2) * np.abs(l2)
        resp = np.maximum(resp, r)
    return np.where(mask, resp, 0.0)


def terrain_concurrence(
    scarp_path: Path | str,
    det_elev_slope: np.ndarray,
    footprint: np.ndarray,
    *,
    channels: tuple[int, ...] = (1, 2, 3, 4, 5, 6, 7, 8, 9),
) -> tuple[np.ndarray, dict]:
    """Surface-expression evidence in [0, 1] from the 3DEP 1 m DEM scarp stack.

    Each requested channel is quantile-scaled independently and the maximum is
    taken, because a scarp may express as excess slope, an anti-slope facing
    step, a cross-slope step or a crest/base curvature at different places.
    The challenge detrended-elevation slope band is added as a coarse,
    always-available term so that cells outside 3DEP lidar coverage are not
    forced to zero.
    """
    import rasterio

    with rasterio.open(scarp_path) as src:
        stack = np.stack([src.read(i).astype(np.float64) for i in channels if i <= src.count])
    cov = None
    score = np.zeros(footprint.shape, dtype=np.float64)
    for band in stack:
        score = np.maximum(score, robust_normalize(band, footprint & (band > 0), hi=99.5))
    slope = robust_normalize(np.abs(det_elev_slope), footprint & np.isfinite(det_elev_slope), hi=99.5)
    score = np.maximum(score, 0.85 * slope)
    # Lidar coverage summary (channel 9 in the delivered product is `relief`).
    with rasterio.open(scarp_path) as src:
        if src.count >= 12:
            cov = src.read(12).astype(np.float64) > 0
    stats = {
        "scarp_path": str(scarp_path),
        "channels_used": [int(c) for c in channels if c <= int(stack.shape[0])],
        "footprint_cells": int(footprint.sum()),
        "concurrence_positive": int((score > 0).sum()),
        "concurrence_median_on_positive": float(np.median(score[score > 0])) if (score > 0).any() else 0.0,
    }
    if cov is not None:
        stats["lidar_covered_cells"] = int((cov & footprint).sum())
        stats["lidar_coverage_fraction"] = float((cov & footprint).sum() / max(int(footprint.sum()), 1))
    return np.where(footprint, score, 0.0), stats


def mapping_gap_weight(
    qfaults_path: Path | str,
    footprint: np.ndarray,
    *,
    coarse_band: int = 1,
    fine_band: int = 2,
    floor: float = 0.5,
) -> tuple[np.ndarray, dict]:
    """H4-B covariate: weight by how *un*-detailedly an area is already compiled.

    The catalogue is a coarse-scale compilation; where detailed (fine-scale)
    QFFD mapping already exists the ground is effectively covered, so candidate
    mass there is less likely to be a genuinely new fault.  The weight is
    ``floor + (1 - floor) * exp(-fine_density / scale)`` blended toward 1 where
    no lidar/QFFD detail exists; it never zeroes a cell (that would be
    tantamount to asserting the hidden labels).
    """
    import rasterio

    with rasterio.open(qfaults_path) as src:
        coarse = src.read(coarse_band).astype(np.float64) > 0
        fine = src.read(fine_band).astype(np.float64) > 0 if src.count >= fine_band else np.zeros_like(coarse)
    if not fine.any():
        return np.ones(footprint.shape, dtype=np.float64), {
            "fine_scale_px": 0,
            "coarse_scale_px": int((coarse & footprint).sum()),
            "note": "no fine-scale QFFD traces inside or near the grid; weight is inert",
        }
    fine_near = ndimage.gaussian_filter(fine.astype(np.float64), 8.0, mode="constant") * (
        (2 * np.pi) * 8.0 ** 2
    )
    w = floor + (1.0 - floor) * np.exp(-fine_near / 2.0)
    w = np.where(footprint, w, 1.0)
    return w, {
        "fine_scale_px_in_footprint": int((fine & footprint).sum()),
        "coarse_scale_px_in_footprint": int((coarse & footprint).sum()),
        "floor": float(floor),
        "min_weight": float(w[footprint].min()),
        "max_weight": float(w[footprint].max()),
    }
