#!/usr/bin/env python3
"""Turn the H45 Euler solution cloud into a continuous, depth-weighted density raster.

The specification being implemented is:

    "Convert that cloud into a continuous raster via kernel-density estimation of solution
     density per pixel, weighted so tight clusters of shallow, mutually-consistent solutions
     score higher than scattered or deep ones, since a real near-surface fault produces the
     former and noise produces the latter."

Weighting
---------
Each accepted Euler solution ``i`` carries a weight

    w_i = w_depth(z_i) * w_fit(r_i) * w_precision(se_i) * w_tight(i) * w_cross(i)

* ``w_depth = exp(-z / D0)``            shallow solutions dominate; ``D0`` is preregistered.
* ``w_fit   = exp(-(r/r0)^2 / 2)``      relative least-squares misfit of the Euler equation.
* ``w_precision = exp(-(se/(f*z))^2 / 2)``   conditional depth standard error.
* ``w_tight = exp(-(mad_i / tol)^2 / 2)``    local depth scatter of neighbouring solutions
  (robust MAD), i.e. *mutual consistency* in depth.
* ``w_cross = (1 + k * n_other_scales) / (1 + k * n_scales)``  a solution corroborated by
  other window sizes / structural indices of the same field is a coherent structure rather
  than a single-window artefact.

Clustering is deliberately *soft*: instead of a hard "keep if it has >= 2 neighbours" rule
(the H4 solver's rule, which deleted isolated but possibly real solutions), the neighbour
count and depth scatter enter as continuous factors.

Rasterisation
-------------
Weights are splatted bilinearly at sub-pixel solution positions and then smoothed with a
separable Gaussian kernel (truncated at 4 sigma), which is a kernel-density estimate on the
100 m lattice.  Everything outside the competition footprint is left as NaN.

Optimisation note: the KDE is a convolution, so it is done with ``scipy.ndimage`` separated
1-D kernels rather than an explicit N^2 kernel matrix.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

from .euler_h45 import SolutionCloud


@dataclass(frozen=True)
class ClusterConfig:
    depth_scale_m: float = 900.0        # D0 for the shallowness weight
    residual_scale: float = 0.30        # r0 for the misfit weight
    depth_se_frac: float = 0.35         # f for the precision weight
    neighbour_radius_m: float = 700.0   # horizontal radius for the consistency statistic
    scatter_tol_m: float = 400.0        # tolerance for the robust local depth scatter
    min_neighbours: int = 2             # neighbours (excluding self) needed for tightness > 0
    cross_scale_bonus: float = 0.75     # k in the cross-scale corroboration factor
    kde_sigma_px: float = 2.0
    kde_truncate: float = 4.0
    saturation: str = "log1p"           # "log1p" or "none": per-family dynamic-range control


def _robust_local_scatter(xy: np.ndarray, depth: np.ndarray, radius: float,
                          k: int = 8) -> tuple[np.ndarray, np.ndarray]:
    """MAD of neighbour depths and neighbour count (excluding self), for every solution."""
    tree = cKDTree(xy)
    counts = np.fromiter((len(p) for p in tree.query_ball_point(xy, r=radius)),
                         dtype=np.int64, count=len(xy)) - 1
    # k nearest neighbours (self included) -> MAD over the other k-1
    kk = max(2, min(k, len(xy)))
    _, idx = tree.query(xy, k=kk, workers=1)
    idx = np.atleast_2d(idx)
    nb = depth[idx]                       # (n, kk)
    med = np.median(nb, axis=1, keepdims=True)
    mad = np.median(np.abs(nb - med), axis=1) * 1.4826
    return mad, counts


def _cross_scale_counts(xy: np.ndarray, tags: np.ndarray, radius: float) -> np.ndarray:
    """How many *distinct other* (window, SI) groups have a solution within ``radius``."""
    tree = cKDTree(xy)
    n = len(xy)
    out = np.zeros(n, dtype=np.int64)
    groups = np.unique(tags)
    for g in groups:
        other = np.flatnonzero(tags != g)
        if other.size == 0:
            continue
        sub = cKDTree(xy[other])
        hits = np.array(len(p) for p in sub.query_ball_point(xy[tags == g], r=radius))
        out[tags == g] = (hits > 0).astype(np.int64)   # at least one neighbour elsewhere
    return out


def cluster_weights(clouds: list[SolutionCloud], cfg: ClusterConfig,
                    dx: float = 100.0) -> list[SolutionCloud]:
    """Attach a ``weight`` array to each cloud (in place) and return them."""
    for cloud in clouds:
        n = len(cloud)
        if n == 0:
            cloud.meta["weight"] = np.zeros(0)
            continue
        z = np.maximum(cloud.depth, 0.0)
        w_depth = np.exp(-z / cfg.depth_scale_m)
        w_fit = np.exp(-0.5 * (cloud.rel_residual / cfg.residual_scale) ** 2)
        denom = np.maximum(cfg.depth_se_frac * np.maximum(z, 50.0), 25.0)
        w_prec = np.exp(-0.5 * (cloud.depth_se / denom) ** 2)
        xy = np.column_stack([cloud.col * dx, cloud.row * dx])
        mad, counts = _robust_local_scatter(xy, cloud.depth, cfg.neighbour_radius_m)
        w_tight = np.exp(-0.5 * (mad / cfg.scatter_tol_m) ** 2)
        w_tight = np.where(counts >= cfg.min_neighbours, w_tight, 0.35 * w_tight)
        w = w_depth * w_fit * w_prec * w_tight
        cloud.meta["weight"] = w
        cloud.meta["neighbour_counts"] = counts
        cloud.meta["local_mad"] = mad
        cloud.meta["xy"] = xy
    return clouds


def cross_family_weight(clouds: list[SolutionCloud], cfg: ClusterConfig,
                        dx: float = 100.0) -> None:
    """Add the cross-scale / cross-index corroboration factor inside each field family."""
    by_field: dict[str, list[int]] = {}
    for i, c in enumerate(clouds):
        by_field.setdefault(c.family, []).append(i)
    for field, idxs in by_field.items():
        idxs = [i for i in idxs if len(clouds[i]) > 0]
        if not idxs:
            continue
        if len(idxs) < 2:
            for i in idxs:
                clouds[i].meta["cross"] = np.ones(len(clouds[i]))
            continue
        xy = np.concatenate([clouds[i].meta["xy"] for i in idxs])
        tags = np.concatenate([np.full(len(clouds[i]), k) for k, i in enumerate(idxs)])
        # number of *other* groups that have at least one solution nearby
        radius = cfg.neighbour_radius_m
        tree = cKDTree(xy)
        groups = np.unique(tags)
        per_group_presence = np.zeros((len(xy), len(groups)), dtype=bool)
        for k, g in enumerate(groups):
            other = np.flatnonzero(tags != g)
            if other.size == 0:
                continue
            sub = cKDTree(xy[other])
            pts = xy[tags == g]
            hits = np.fromiter((len(p) for p in sub.query_ball_point(pts, r=radius)),
                               dtype=np.int64, count=len(pts))
            per_group_presence[tags == g, k] = hits > 0
        present = per_group_presence.sum(axis=1) - 1      # exclude own group
        present = np.maximum(present, 0)
        n_other = max(len(groups) - 1, 1)
        factor = (1.0 + cfg.cross_scale_bonus * present) / (1.0 + cfg.cross_scale_bonus * n_other)
        offset = 0
        for i in idxs:
            n_i = len(clouds[i])
            clouds[i].meta["cross"] = factor[offset:offset + n_i]
            clouds[i].meta["weight"] = clouds[i].meta["weight"] * clouds[i].meta["cross"]
            offset += n_i


def splat(cloud: SolutionCloud, shape: tuple[int, int], dx: float = 100.0) -> np.ndarray:
    """Bilinear mass-conserving splat of the solution weights onto the 100 m lattice."""
    acc = np.zeros(shape, dtype=np.float64)
    if len(cloud) == 0:
        return acc
    w = cloud.meta.get("weight", np.ones(len(cloud)))
    r = np.clip(cloud.row, 0.0, shape[0] - 1.001)
    c = np.clip(cloud.col, 0.0, shape[1] - 1.001)
    r0 = np.floor(r).astype(np.int64)
    c0 = np.floor(c).astype(np.int64)
    fr = r - r0
    fc = c - c0
    r1 = np.minimum(r0 + 1, shape[0] - 1)
    c1 = np.minimum(c0 + 1, shape[1] - 1)
    np.add.at(acc, (r0, c0), w * (1 - fr) * (1 - fc))
    np.add.at(acc, (r0, c1), w * (1 - fr) * fc)
    np.add.at(acc, (r1, c0), w * fr * (1 - fc))
    np.add.at(acc, (r1, c1), w * fr * fc)
    return acc


def kde(field: np.ndarray, sigma_px: float = 2.0, truncate: float = 4.0) -> np.ndarray:
    """Separable Gaussian kernel-density estimate on the lattice (constant-zero padding)."""
    return ndimage.gaussian_filter(field, sigma=sigma_px, truncate=truncate,
                                   mode="constant", cval=0.0)


def normalise(field: np.ndarray, valid: np.ndarray, mode: str = "max",
              percentile: float = 99.9) -> np.ndarray:
    """Scale a non-negative field into [0, 1]."""
    out = np.where(valid, field, 0.0)
    out = np.maximum(out, 0.0)
    if mode == "max":
        m = out.max()
    elif mode == "percentile":
        v = out[out > 0]
        m = float(np.percentile(v, percentile)) if v.size else 0.0
    else:
        raise ValueError(mode)
    if m <= 0:
        return np.zeros_like(out)
    return np.clip(out / m, 0.0, 1.0)
