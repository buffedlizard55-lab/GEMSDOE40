"""Measurement layer for the organizer's distance-weighted Tversky index.

This is the fast, auditable route to the exact components of the metric on a
12.3-million-cell grid.  It is deliberately redundant with
:mod:`gemsdoe40.research_metric` (which is pinned to the organizer's worked
example) so that two independent implementations must agree: the tests check
them against each other on random fields.

Definitions (problem description, page 967), alpha = 0.2, beta = 0.8, R = 3 px:

    TPw = sum_{g in G} max_{x: d(x,g) <= R} p(x) k(d(x,g))
    FPw = sum_x p(x) [1 - max_{g in G} k(d(x,g))]
    FNw = sum_g [1 - max_{x: d(x,g) <= R} p(x) k(d(x,g))]  =  N - TPw
    DTI = TPw / (TPw + 0.2 FPw + 0.8 FNw + eps)

The fast route uses the identity

    TPw = sum_{g in G} K[g],   K[y, x] = max_{|dy|,|dx| <= R, d <= R} p[y+dy, x+dx] k(d)

so one (2R+1)^2 shift-max pass gives TPw for any truth mask, and one Euclidean
distance transform of the truth mask gives FPw for any prediction.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

ALPHA = 0.2
BETA = 0.8
RADIUS_PX = 3.0


def _offsets(radius_px: float) -> list[tuple[int, int, float]]:
    rad = int(np.floor(radius_px))
    out: list[tuple[int, int, float]] = []
    for dy in range(-rad, rad + 1):
        for dx in range(-rad, rad + 1):
            d = float(np.hypot(dx, dy))
            if d <= radius_px:
                out.append((dy, dx, 1.0 - d / radius_px))
    return out


def kernel_field(prediction: np.ndarray, radius_px: float = RADIUS_PX) -> np.ndarray:
    """``K[y, x] = max over the R-neighbourhood of p * triangular kernel``.

    Float32, same shape as the prediction.  ``K[g]`` is exactly the credit that
    ground-truth pixel ``g`` receives from this prediction, so
    ``TPw = K[truth].sum()``.
    """
    p = np.asarray(prediction, dtype=np.float32)
    if p.ndim != 2:
        raise ValueError("prediction must be 2-D")
    h, w = p.shape
    out = np.zeros((h, w), dtype=np.float32)
    for dy, dx, wt in _offsets(radius_px):
        shifted = np.zeros_like(p)
        ys = slice(max(0, dy), min(h, h + dy))
        yd = slice(max(0, -dy), min(h, h - dy))
        xs = slice(max(0, dx), min(w, w + dx))
        xd = slice(max(0, -dx), min(w, w - dx))
        shifted[ys, xs] = p[yd, xd]
        np.maximum(out, shifted * np.float32(wt), out=out)
    return out


def kernel_from_truth(truth: np.ndarray, radius_px: float = RADIUS_PX) -> np.ndarray:
    """``max_g k(d(x, g))`` per pixel: the factor that scales a prediction into FPW."""
    g = np.asarray(truth, dtype=bool)
    if not g.any():
        return np.zeros(g.shape, dtype=np.float32)
    d = ndimage.distance_transform_edt(~g)
    return np.maximum(1.0 - d / float(radius_px), 0.0).astype(np.float32)


def credit_components(
    prediction: np.ndarray,
    truth: np.ndarray,
    *,
    valid: np.ndarray | None = None,
    mask: np.ndarray | None = None,
    radius_px: float = RADIUS_PX,
    kernel_field_cache: np.ndarray | None = None,
    kernel_truth_cache: np.ndarray | None = None,
) -> dict:
    """Exact TPw / FPw / FNw / N / mass / DTI for one raster.

    ``valid`` restricts the scored domain (default: everywhere).
    ``mask`` is the pixel-exact known-fault mask: masked cells are excluded from
    both the prediction and the truth, exactly as the organizer confirmed.
    """
    p = np.asarray(prediction, dtype=np.float32)
    g = np.asarray(truth, dtype=bool)
    if p.shape != g.shape:
        raise ValueError("prediction and truth must share a shape")
    score_all = np.ones(p.shape, dtype=bool) if valid is None else np.asarray(valid, dtype=bool)
    if mask is not None:
        m = np.asarray(mask, dtype=bool)
        p = np.where(m, 0.0, p).astype(np.float32)
        g = g & ~m
    g = g & score_all
    p = np.where(score_all, p, 0.0).astype(np.float32)

    K = kernel_field(p, radius_px) if kernel_field_cache is None else kernel_field_cache
    if mask is not None:
        K = np.where(np.asarray(mask, dtype=bool), 0.0, K)
    tp = float(K[g].sum(dtype=np.float64)) if g.any() else 0.0
    kt = kernel_from_truth(g, radius_px) if kernel_truth_cache is None else kernel_truth_cache
    if mask is not None:
        kt = np.where(np.asarray(mask, dtype=bool), 0.0, kt)
    pred_credit = float((p * kt).sum(dtype=np.float64))
    mass = float(p.sum(dtype=np.float64))
    fp = mass - pred_credit
    n_truth = int(g.sum())
    fn = n_truth - tp
    denom = tp + ALPHA * fp + BETA * fn
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "pred_credit": pred_credit,
        "mass": mass,
        "n_truth": n_truth,
        "dti": tp / denom if denom > 0 else 0.0,
        "credit_per_mass": tp / mass if mass > 0 else 0.0,
    }


def blocked_components(
    prediction: np.ndarray,
    truth: np.ndarray,
    *,
    valid: np.ndarray,
    mask: np.ndarray | None = None,
    n_rows: int = 4,
    n_cols: int = 6,
    guard: int = 3,
    radius_px: float = RADIUS_PX,
) -> dict:
    """Pooled + per-block DTI on disjoint block interiors with a guard band.

    The guard removes cross-block credit transfer, so a block's score depends
    only on predictions and truth inside that block.  Mirrors
    :func:`gemsdoe40.research_metric.spatial_block_components` while using the
    fast kernels above.
    """
    p = np.asarray(prediction, dtype=np.float32)
    g = np.asarray(truth, dtype=bool)
    v = np.asarray(valid, dtype=bool)
    h, w = p.shape
    y_edges = np.linspace(0, h, n_rows + 1, dtype=int)
    x_edges = np.linspace(0, w, n_cols + 1, dtype=int)
    scored = np.zeros((h, w), dtype=bool)
    boxes: list[tuple[int, int, slice, slice]] = []
    for iy in range(n_rows):
        r0, r1 = int(y_edges[iy]), int(y_edges[iy + 1])
        rr0 = r0 + (guard if iy > 0 else 0)
        rr1 = r1 - (guard if iy < n_rows - 1 else 0)
        for ix in range(n_cols):
            c0, c1 = int(x_edges[ix]), int(x_edges[ix + 1])
            cc0 = c0 + (guard if ix > 0 else 0)
            cc1 = c1 - (guard if ix < n_cols - 1 else 0)
            rs, cs = slice(rr0, rr1), slice(cc0, cc1)
            scored[rs, cs] = True
            boxes.append((iy, ix, rs, cs))
    pooled = credit_components(p, g, valid=scored & v, mask=mask, radius_px=radius_px)
    blocks = []
    for iy, ix, rs, cs in boxes:
        sub = np.zeros((h, w), dtype=bool)
        sub[rs, cs] = True
        blk = credit_components(p, g, valid=sub & v, mask=mask, radius_px=radius_px)
        blk.update({"row": iy, "col": ix})
        blocks.append(blk)
    scores = [b["dti"] for b in blocks]
    return {
        "pooled": pooled,
        "blocks": blocks,
        "mean_block_dti": float(np.mean(scores)) if scores else 0.0,
        "folds_positive_vs_zero": int(sum(s > 0 for s in scores)),
    }


def same_mass_random(
    prediction: np.ndarray,
    valid: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Binary control with the same emitted mass, drawn uniformly from ``valid``."""
    n = int(np.count_nonzero(np.asarray(prediction) > 0))
    idx = np.flatnonzero(np.asarray(valid, dtype=bool).ravel())
    n = min(n, idx.size)
    out = np.zeros(prediction.shape, dtype=np.float32)
    if n:
        out.ravel()[rng.choice(idx, size=n, replace=False)] = 1.0
    return out
