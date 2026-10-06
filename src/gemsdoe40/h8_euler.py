"""H8: lineament-weighted, cross-family-corroborated Euler depth-clustering.

Construction (all of it from ``training_features.tif``; no prior prediction and
no holdout raster is read):

1. SI = 0 contact Euler deconvolution (Reid et al. 1990, eq. 2 with the
   arbitrary offset A) on multi-window patches of
     * magnetics: reduced-to-pole total field (band ``rtp``);
     * gravity:   vertical derivative of the isostatic gravity anomaly
                  (band ``iso_grav_anom``, differentiated once in the FFT),
   which is the pairing for which the fitted field is homogeneous of degree 0
   at a two-dimensional contact edge.  Reid & Thurston (2014) show a *finite*
   gravity step has SI = -1 and needs a generalised treatment, so the gravity
   arm is used as a local top-edge approximation, exactly as H7/H4 recorded.
2. QC by design condition, relative residual, depth range and depth standard
   error (the frozen H4 gates, re-parameterised here -- the H4 constants in
   :mod:`gemsdoe40.contact_euler` are NOT modified, so H4 stays reproducible).
3. Weighting, in the order the brief names them:
     * shallow-first  exp(-depth / depth_decay);
     * tight clusters: cross-window depth consensus within a 300 m xy radius;
     * mutual consistency: local *lineament* coherence from the anisotropy of
       the solution cloud (a fault is a line, not a blob);
     * independent corroboration: a solution is boosted only when the other
       family independently places a solution at the same place and depth.
4. Kernel-density estimation of the weighted cloud per pixel, normalised to
   [0, 1] by its own high percentile.
5. Metric-optimal emission: greedy value-ranked max-separation thinning at
   2.8 px (just inside the metric's own 3 px kernel radius) under a fixed mass
   budget, as a *continuous* field normalised to a maximum of 1.0.

Nothing here is a score.  Every number this module produces is a measurement on
a named surrogate or a model, and is labelled as such by the runner.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

# --------------------------------------------------------------------------- #
# configuration
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class FamilyConfig:
    """One potential-field family: a band, its derivative order, windows."""

    name: str
    continuation_m: float
    derivative_order: int
    windows: tuple[int, ...]
    edge_guard_m: float
    max_depth_m: float


FAMILIES = (
    # magnetic reduced-to-pole total field: SI = 0 is the contact index exactly
    FamilyConfig("rtp", 150.0, 0, (7, 11, 15, 21), 700.0, 3000.0),
    # isostatic gravity, differentiated once: the contact-like arm (approximation)
    FamilyConfig("iso_grav_anom", 500.0, 1, (11, 15, 21, 31), 1500.0, 5000.0),
)

SETTINGS: dict[str, Any] = {
    "si": 0,
    "fit_contact_offset": True,
    "resolution_m": 100.0,
    "fft_pad": 128,
    "stride": 4,
    "max_design_condition": 200.0,
    "max_relative_residual": 0.25,
    "min_effective_depth_m": 50.0,
    "max_depth_se_fraction": 0.20,
    "depth_se_floor_m": 75.0,
    # cluster / coherence / corroboration radii (metres)
    "cluster_xy_m": 300.0,
    "cluster_depth_floor_m": 200.0,
    "cluster_depth_fraction": 0.35,
    "min_other_neighbours": 3,
    "coherence_xy_m": 500.0,
    "min_coherence_neighbours": 5,
    "corroboration_xy_m": 400.0,
    "corroboration_depth_m": 600.0,
    "uncorroborated_factor": 0.5,
    "depth_decay_m": 1200.0,
    # density field and emission
    "kde_sigma_px": 1.5,
    "kde_truncate": 4.0,
    "family_percentile": 99.0,
    "emission_spacing_px": 2.8,
    "mass_budget": 45000,
}

CLOUD_DTYPE = np.dtype([
    ("row", "f8"), ("col", "f8"), ("depth_m", "f8"), ("depth_se_m", "f8"),
    ("residual", "f8"), ("condition", "f8"), ("offset", "f8"),
    ("rank", "u1"), ("window", "u1"), ("family_id", "u1"), ("weight", "f8"),
    ("neighbours", "i4"), ("depth_mad_m", "f8"), ("coherence", "f8"),
    ("corroborated", "u1"),
])


def config_dict() -> dict[str, Any]:
    return {"settings": SETTINGS, "families": [asdict(f) for f in FAMILIES]}


# --------------------------------------------------------------------------- #
# weighting
# --------------------------------------------------------------------------- #


def depth_consensus_weights(cloud: np.ndarray, *, xy_radius_m: float,
                            depth_floor_m: float, depth_fraction: float,
                            min_neighbours: int, depth_decay_m: float,
                            cell_m: float = 100.0) -> np.ndarray:
    """Shallow-first weight gated on cross-window depth consensus.

    A solution keeps a non-zero weight only when at least ``min_neighbours``
    other solutions lie within ``xy_radius_m`` **and** within a depth tolerance,
    and at least one of them comes from a *different* window size.  The
    magnitude combines shallowness, fit residual and depth standard error.
    """
    result = cloud.copy()
    if not len(cloud):
        return result
    tree = cKDTree(np.column_stack((cloud["row"], cloud["col"])) * cell_m)
    for start in range(0, len(cloud), 512):
        points = tree.data[start:start + 512]
        neighbours = tree.query_ball_point(points, xy_radius_m)
        for offset, ids in enumerate(neighbours):
            i = start + offset
            ids = np.asarray(ids, dtype=int)
            tolerance = max(depth_floor_m, depth_fraction * cloud["depth_m"][i])
            ids = ids[(ids != i) & (np.abs(cloud["depth_m"][ids] - cloud["depth_m"][i]) <= tolerance)]
            n = len(ids)
            result["neighbours"][i] = n
            if n < min_neighbours or not np.any(cloud["window"][ids] != cloud["window"][i]):
                continue
            median = np.median(cloud["depth_m"][ids])
            mad = 1.4826 * np.median(np.abs(cloud["depth_m"][ids] - median))
            depth = cloud["depth_m"][i]
            residual = cloud["residual"][i]
            se = cloud["depth_se_m"][i]
            weight = (np.exp(-depth / depth_decay_m)
                      * np.exp(-0.5 * (residual / SETTINGS["max_relative_residual"]) ** 2)
                      * np.exp(-0.5 * (se / max(depth, 100.0) / SETTINGS["max_depth_se_fraction"]) ** 2)
                      * (n / (n + 3.0))
                      * np.exp(-0.5 * (mad / tolerance) ** 2))
            result["weight"][i] = weight
            result["depth_mad_m"][i] = mad
    return result


def lineament_coherence(cloud: np.ndarray, *, xy_radius_m: float,
                        min_neighbours: int, cell_m: float = 100.0) -> np.ndarray:
    """Local anisotropy of the solution cloud: 1 on a line, 0 for a blob.

    For every solution with non-zero weight we take the neighbouring solutions
    inside ``xy_radius_m`` and compute the 2x2 covariance of their offsets.  A
    fault trace produces solutions strung out along a line, so the smaller
    eigenvalue is small relative to the larger one; a noise field produces
    isotropic scatter.  ``coherence = 1 - lambda_2 / lambda_1``.
    """
    coherence = np.zeros(len(cloud), dtype=np.float64)
    active = np.flatnonzero(cloud["weight"] > 0)
    if not len(active):
        return coherence
    tree = cKDTree(np.column_stack((cloud["row"], cloud["col"])) * cell_m)
    for start in range(0, len(active), 512):
        chunk = active[start:start + 512]
        points = np.column_stack((cloud["row"][chunk], cloud["col"][chunk])) * cell_m
        neighbours = tree.query_ball_point(points, xy_radius_m)
        for offset, ids in enumerate(neighbours):
            i = chunk[offset]
            ids = np.asarray(ids, dtype=int)
            if len(ids) < min_neighbours:
                continue
            delta = np.column_stack((cloud["row"][ids], cloud["col"][ids])) * cell_m - points[offset]
            cov = (delta.T @ delta) / len(ids)
            eigenvalues = np.linalg.eigvalsh(cov)
            if eigenvalues[1] <= 0:
                continue
            coherence[i] = float(np.clip(1.0 - eigenvalues[0] / eigenvalues[1], 0.0, 1.0))
    return coherence


def cross_family_corroboration(magnetic: np.ndarray, gravity: np.ndarray, *,
                               xy_radius_m: float, depth_m: float,
                               cell_m: float = 100.0) -> tuple[np.ndarray, np.ndarray]:
    """1 where the other family independently places a solution nearby in x, y, z."""
    flags = [np.zeros(len(c), dtype=np.uint8) for c in (magnetic, gravity)]
    trees = [cKDTree(np.column_stack((c["row"], c["col"])) * cell_m) if len(c) else None
             for c in (magnetic, gravity)]
    for k, (cloud, other, tree) in enumerate(((magnetic, gravity, trees[1]), (gravity, magnetic, trees[0]))):
        if tree is None or not len(cloud):
            continue
        for start in range(0, len(cloud), 512):
            chunk = slice(start, min(start + 512, len(cloud)))
            points = np.column_stack((cloud["row"][chunk], cloud["col"][chunk])) * cell_m
            neigh = tree.query_ball_point(points, xy_radius_m)
            for offset, ids in enumerate(neigh):
                i = start + offset
                ids = np.asarray(ids, dtype=int)
                if not len(ids):
                    continue
                if np.any(np.abs(other["depth_m"][ids] - cloud["depth_m"][i]) <= depth_m):
                    flags[k][i] = 1
    return flags[0], flags[1]


def combine_weights(cloud: np.ndarray) -> np.ndarray:
    """Fold coherence and corroboration into the depth-consensus weight."""
    coherence = np.sqrt(np.clip(cloud["coherence"], 0.0, 1.0))
    factor = np.where(cloud["corroborated"] > 0, 1.0, SETTINGS["uncorroborated_factor"])
    weight = np.where(cloud["weight"] > 0, cloud["weight"] * coherence * factor, 0.0)
    if not np.isfinite(weight).all():
        raise ValueError("non-finite combined weight")
    return weight.astype(np.float64)


# --------------------------------------------------------------------------- #
# density field and emission
# --------------------------------------------------------------------------- #


def splat(rows: np.ndarray, cols: np.ndarray, weights: np.ndarray,
          shape: tuple[int, int]) -> np.ndarray:
    """Bilinearly deposit point mass at sub-pixel positions."""
    result = np.zeros(shape, dtype=np.float64)
    r0, c0 = np.floor(rows).astype(int), np.floor(cols).astype(int)
    fr, fc = rows - r0, cols - c0
    for dr, rw in ((0, 1 - fr), (1, fr)):
        for dc, cw in ((0, 1 - fc), (1, fc)):
            rr, cc = r0 + dr, c0 + dc
            good = (rr >= 0) & (rr < shape[0]) & (cc >= 0) & (cc < shape[1])
            np.add.at(result, (rr[good], cc[good]), (weights * rw * cw)[good])
    return result


def kde_field(cloud: np.ndarray, weight: np.ndarray, footprint: np.ndarray) -> tuple[np.ndarray, dict]:
    """Weighted per-pixel kernel density, normalised to [0, 1] inside the footprint."""
    retained = weight > 0
    if not retained.any():
        raise ValueError("no retained solutions; refusing to backfill")
    raw = splat(cloud["row"][retained], cloud["col"][retained], weight[retained], footprint.shape)
    kde = ndimage.gaussian_filter(raw, sigma=SETTINGS["kde_sigma_px"],
                                  truncate=SETTINGS["kde_truncate"], mode="constant", cval=0.0)
    kde[~footprint] = 0.0
    support = kde[footprint & (kde > 0)]
    scale = float(np.percentile(support, SETTINGS["family_percentile"]))
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("degenerate KDE normaliser")
    field = np.clip(kde / scale, 0.0, 1.0).astype(np.float32)
    stats = {
        "solutions": int(len(cloud)),
        "clustered_solutions": int(retained.sum()),
        "positive_cells": int(len(support)),
        "normaliser": scale,
        "weight_sum": float(weight.sum()),
        "depth_m_p05_p50_p95": np.percentile(cloud["depth_m"][retained], [5, 50, 95]).tolist(),
        "coherence_p50": float(np.median(cloud["coherence"][retained])),
        "corroborated_fraction": float((cloud["corroborated"][retained] > 0).mean()),
    }
    return field, stats
