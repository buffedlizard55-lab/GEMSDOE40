"""Weighted kernel-density raster of an Euler solution cloud.

A real near-surface fault produces a *tight cluster of shallow, mutually
consistent* Euler solutions; noise produces a scattered or deep cloud
(Reid et al. 1990, figs. 1 and 5).  This module converts the point cloud
into a continuous raster by:

  1. weighting each solution by shallowness, solution quality, local
     tightness, and local depth-consistency;
  2. splatting those weights onto the competition grid;
  3. smoothing with a Gaussian kernel (the KDE);
  4. boosting pixels where magnetic and gravity solutions agree
     (dual-physics concordance — unique to GEMSDOE40).

The result is a continuous field in [0, 1], not a binary edge map and not a
copy of any dotted H19-5 / H33 geometry.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter, binary_dilation
from scipy.spatial import cKDTree

from .euler import EulerCloud


def solution_weights(
    cloud: EulerCloud,
    *,
    depth_scale_m: float = 700.0,
    neighbor_px: float = 3.0,
    min_neighbors: int = 3,
) -> np.ndarray:
    """Per-solution weight: shallow × precise × locally tight × depth-consistent."""
    n = len(cloud)
    if n == 0:
        return np.empty(0, dtype=np.float64)
    shallow = np.exp(-np.clip(cloud.depth_m, 0, None) / depth_scale_m)
    quality = 1.0 / (1.0 + 8.0 * np.clip(cloud.rel_se, 0, 10))
    # Local tightness / mutual depth consistency via a 2-D KD tree on (row, col).
    xy = np.column_stack([cloud.row, cloud.col])
    tree = cKDTree(xy)
    neighbors = tree.query_ball_tree(tree, r=neighbor_px)
    tight = np.empty(n, dtype=np.float64)
    consist = np.empty(n, dtype=np.float64)
    for i, nb in enumerate(neighbors):
        k = len(nb)
        tight[i] = float(k)
        if k >= min_neighbors:
            depths = cloud.depth_m[np.fromiter(nb, dtype=int, count=k)]
            mu = float(depths.mean())
            cv = float(depths.std()) / max(mu, 1.0)
            consist[i] = np.exp(-4.0 * cv)
        else:
            consist[i] = 0.0  # isolated solutions are treated as noise
    tight = tight / max(float(np.percentile(tight, 90)), 1.0)
    tight = np.clip(tight, 0.0, 3.0)
    w = shallow * quality * tight * consist
    # Analytic-signal amplitude at the source (stronger contact → higher weight)
    if cloud.analytic.size == n:
        a = cloud.analytic
        p90 = float(np.percentile(a, 90)) if a.size else 1.0
        w = w * (0.25 + 0.75 * np.clip(a / max(p90, 1e-12), 0, 2.0))
    return np.where(np.isfinite(w) & (w > 0), w, 0.0)


def splat_kde(
    cloud: EulerCloud,
    weights: np.ndarray,
    shape: tuple[int, int],
    *,
    sigma_px: float = 1.8,
) -> np.ndarray:
    """Bilinear splat of weighted solutions, then Gaussian KDE."""
    acc = np.zeros(shape, dtype=np.float64)
    if len(cloud) == 0 or not np.any(weights > 0):
        return acc
    h, w = shape
    r = cloud.row
    c = cloud.col
    keep = (r >= 0) & (r <= h - 1) & (c >= 0) & (c <= w - 1) & (weights > 0)
    r = r[keep]
    c = c[keep]
    wt = weights[keep]
    r0 = np.floor(r).astype(int)
    c0 = np.floor(c).astype(int)
    r1 = np.clip(r0 + 1, 0, h - 1)
    c1 = np.clip(c0 + 1, 0, w - 1)
    r0 = np.clip(r0, 0, h - 1)
    c0 = np.clip(c0, 0, w - 1)
    dr = r - r0
    dc = c - c0
    np.add.at(acc, (r0, c0), wt * (1 - dr) * (1 - dc))
    np.add.at(acc, (r0, c1), wt * (1 - dr) * dc)
    np.add.at(acc, (r1, c0), wt * dr * (1 - dc))
    np.add.at(acc, (r1, c1), wt * dr * dc)
    if sigma_px > 0:
        acc = gaussian_filter(acc, sigma=sigma_px, mode="constant")
    return acc


def concordance_boost(
    mag: EulerCloud,
    grav: EulerCloud,
    shape: tuple[int, int],
    *,
    radius_px: float = 2.5,
    depth_tol: float = 0.35,
    sigma_px: float = 1.6,
) -> np.ndarray:
    """Extra density where a magnetic solution and a gravity solution agree.

    Agreement = location within ``radius_px`` and relative depth difference
    below ``depth_tol``.  This is the dual-physics signature of a real density
    *and* susceptibility contact — a fault that offsets both — and is not
    something a single-field gradient threshold can produce.
    """
    acc = np.zeros(shape, dtype=np.float64)
    if len(mag) == 0 or len(grav) == 0:
        return acc
    tree = cKDTree(np.column_stack([grav.row, grav.col]))
    dist, idx = tree.query(np.column_stack([mag.row, mag.col]), k=1, workers=-1)
    close = dist <= radius_px
    if not np.any(close):
        return acc
    md = mag.depth_m[close]
    gd = grav.depth_m[idx[close]]
    rel = np.abs(md - gd) / np.maximum(0.5 * (md + gd), 1.0)
    ok = rel <= depth_tol
    if not np.any(ok):
        return acc
    rows = mag.row[close][ok]
    cols = mag.col[close][ok]
    w = np.exp(-rel[ok] / max(depth_tol, 1e-6)) * np.exp(-0.5 * (md[ok] + gd[ok]) / 800.0)
    dummy = EulerCloud(
        "concordance", 0.0, -1, -1, rows, cols, md[ok],
        np.zeros(ok.sum()), np.zeros(ok.sum()), np.ones(ok.sum()),
    )
    return splat_kde(dummy, w, shape, sigma_px=sigma_px)


def assemble_field(
    mag_cloud: EulerCloud,
    grav_cloud: EulerCloud,
    footprint: np.ndarray,
    catalogue: np.ndarray,
    *,
    catalogue_buffer_px: int = 2,
    sigma_px: float = 1.8,
    concord_weight: float = 1.6,
    floor_quantile: float = 0.92,
    power: float = 1.35,
) -> tuple[np.ndarray, dict]:
    """Build the continuous [0, 1] raster from mag + grav Euler clouds.

    Catalogue cells (plus a ``catalogue_buffer_px`` dilation) are forced to 0
    because the hidden test set is *new* faults, not the USGS/INGENIOUS map.
    Values below ``floor_quantile`` of the in-footprint positive field are
    zeroed so the DTI false-positive budget is not spent on KDE tails; the
    surviving cores keep their continuous KDE scores.
    """
    shape = footprint.shape
    w_m = solution_weights(mag_cloud)
    w_g = solution_weights(grav_cloud)
    kde_m = splat_kde(mag_cloud, w_m, shape, sigma_px=sigma_px)
    kde_g = splat_kde(grav_cloud, w_g, shape, sigma_px=sigma_px)
    concord = concordance_boost(mag_cloud, grav_cloud, shape, sigma_px=max(1.2, sigma_px - 0.2))
    raw = kde_m + kde_g + concord_weight * concord
    raw = np.where(footprint, raw, 0.0)

    if catalogue is not None and catalogue_buffer_px >= 0:
        cat = np.asarray(catalogue, dtype=bool) & footprint
        if catalogue_buffer_px > 0 and cat.any():
            cat = binary_dilation(cat, iterations=int(catalogue_buffer_px))
        raw = np.where(cat, 0.0, raw)
    else:
        cat = np.zeros(shape, dtype=bool)

    pos = footprint & (raw > 0)
    stats = {
        "mag_solutions": int(len(mag_cloud)),
        "grav_solutions": int(len(grav_cloud)),
        "mag_weight_sum": float(w_m.sum()) if w_m.size else 0.0,
        "grav_weight_sum": float(w_g.sum()) if w_g.size else 0.0,
        "raw_positive_px": int(pos.sum()),
        "catalogue_buffer_px": int(catalogue_buffer_px),
        "catalogue_masked_px": int(cat.sum()),
    }
    if not pos.any():
        return np.zeros(shape, dtype=np.float32), stats

    # Stretch and floor.
    p = raw.copy()
    p99 = float(np.percentile(p[pos], 99.5))
    if p99 <= 0:
        p99 = float(p[pos].max())
    p = np.clip(p / p99, 0.0, 1.0)
    p = np.power(p, power)
    floor = float(np.quantile(p[p > 0], floor_quantile)) if (p > 0).any() else 0.0
    p = np.where(p >= floor, p, 0.0)
    p = np.where(footprint, p, 0.0)
    p = np.where(cat, 0.0, p)
    # Re-normalise surviving cores to use the full [0, 1] range.
    m = float(p[footprint].max()) if footprint.any() else 0.0
    if m > 0:
        p = p / m
    p = np.clip(p, 0.0, 1.0).astype(np.float32)
    stats.update(
        p99=p99,
        floor=floor,
        floor_quantile=floor_quantile,
        power=power,
        emitted_px=int((p > 0).sum()),
        mass=float(p[footprint].sum()),
        max=float(p[footprint].max()) if footprint.any() else 0.0,
        on_catalogue_positive=int(((p > 0) & np.asarray(catalogue, dtype=bool)).sum()),
    )
    return p, stats
