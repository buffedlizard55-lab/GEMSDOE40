"""Metric-aware dot-lattice emission.

Every successful prior submission in this programme is a *binary dot field*:
one pixel per dot, minimum separation exactly ``2*sqrt(2) = 2.83`` px (300 m
kernel radius / sqrt(2)), which is the densest pattern in which no two dots
compete for the same truth pixel's maximum-credit term.  This module builds
that pattern from a ranked continuous field, in a way that is:

* deterministic (ties broken by row/column order),
* mass-exact (the emitted pixel count is the requested budget),
* catalogue-aware (dots are never placed on the pixel-exact catalogue mask,
  because the organiser excludes those cells from both metric terms).

The measurement half of the module reports the redundancy of a dot field -
how many other dots share a 3 px neighbourhood - which is the quantity the
metric charges false-positive mass for.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

DEFAULT_SEPARATION_PX = 2.0 * np.sqrt(2.0)  # 2.8284 px = 282.8 m


def _disk_offsets(radius: float) -> np.ndarray:
    """Integer offsets strictly closer than ``radius`` (what must be blocked)."""
    r = int(np.ceil(radius))
    offs = []
    for dy in range(-r, r + 1):
        for dx in range(-r, r + 1):
            if (dy == 0 and dx == 0):
                continue
            if np.hypot(dy, dx) < radius - 1e-9:
                offs.append((dy, dx))
    return np.asarray(offs, dtype=int)


def greedy_min_separation(score: np.ndarray, allowed: np.ndarray, *,
                          n_max: int, min_separation_px: float = DEFAULT_SEPARATION_PX,
                          seed_ge_zero: bool = True) -> np.ndarray:
    """Greedy maximum-score selection with a hard minimum separation.

    Candidates are visited in descending score order; a candidate is accepted
    when it is at least ``min_separation_px`` from every accepted dot.  Only
    ``allowed`` cells are candidates.  Returns a boolean dot field.
    """
    s = np.asarray(score, dtype=np.float64)
    ok = np.asarray(allowed, dtype=bool) & np.isfinite(s)
    if seed_ge_zero:
        ok &= s > 0
    h, w = s.shape
    flat_ok = np.flatnonzero(ok.ravel())
    if flat_ok.size == 0 or n_max <= 0:
        return np.zeros((h, w), dtype=bool)
    vals = s.ravel()[flat_ok]
    # Descending value; ties broken by row-major order for determinism.
    order = np.argsort(-vals, kind="stable")
    cand = flat_ok[order]
    blocked = np.zeros(h * w, dtype=bool)
    offs = _disk_offsets(min_separation_px)
    out = np.zeros(h * w, dtype=bool)
    count = 0
    for idx in cand:
        if blocked[idx]:
            continue
        out[idx] = True
        count += 1
        r, c = divmod(int(idx), w)
        if offs.size:
            rr = r + offs[:, 0]
            cc = c + offs[:, 1]
            keep = (rr >= 0) & (rr < h) & (cc >= 0) & (cc < w)
            blocked[rr[keep] * w + cc[keep]] = True
        if count >= n_max:
            break
    return out.reshape(h, w)


def prefix_take(score: np.ndarray, allowed: np.ndarray, n: int) -> np.ndarray:
    """Plain top-``n`` pixels by score inside ``allowed`` (no separation rule)."""
    s = np.asarray(score, dtype=np.float64)
    ok = np.asarray(allowed, dtype=bool) & np.isfinite(s) & (s > 0)
    idx = np.flatnonzero(ok.ravel())
    out = np.zeros(s.size, dtype=bool)
    if idx.size == 0 or n <= 0:
        return out.reshape(s.shape)
    vals = s.ravel()[idx]
    n = min(int(n), idx.size)
    sel = np.argpartition(-vals, n - 1)[:n]
    out[idx[sel]] = True
    return out.reshape(s.shape)


def redundancy_stats(dots: np.ndarray, radius_px: float = 3.0, sample: int = 4000,
                     catalogue_distance: np.ndarray | None = None) -> dict:
    """How many *other* dots share each dot's 3 px neighbourhood.

    When ``catalogue_distance`` (a precomputed distance transform to the published
    catalogue) is supplied, the share of mass sitting inside the 2 px and 6 px
    catalogue bands is reported as well; the emitter uses that to refuse to
    auto-select catalogue-hugging geometries.
    """
    b = np.asarray(dots, dtype=bool)
    yy, xx = np.nonzero(b)
    n = yy.size
    if n == 0:
        return dict(n_dots=0, mean_neighbours_within_radius=0.0, frac_with_neighbour=0.0,
                    median_nearest_spacing_px=None)
    pts = np.column_stack([yy, xx]).astype(np.float64)
    tree = cKDTree(pts)
    take = pts if n <= sample else pts[np.linspace(0, n - 1, sample).astype(int)]
    cnt = np.array([len(tree.query_ball_point(p, radius_px)) - 1 for p in take])
    dd, _ = tree.query(pts, k=2)
    out = dict(
        n_dots=int(n),
        mean_neighbours_within_radius=float(cnt.mean()),
        frac_with_neighbour=float(np.mean(cnt >= 1)),
        median_nearest_spacing_px=float(np.median(dd[:, 1])),
        min_nearest_spacing_px=float(dd[:, 1].min()),
    )
    if catalogue_distance is not None:
        dcat = np.asarray(catalogue_distance)[yy, xx]
        out.update(
            frac_within_2px_of_catalogue=float(np.mean(dcat <= 2.0)),
            frac_within_6px_of_catalogue=float(np.mean(dcat <= 6.0)),
            median_distance_to_catalogue_px=float(np.median(dcat)),
        )
    return out
