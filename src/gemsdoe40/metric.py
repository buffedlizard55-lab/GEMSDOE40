"""Transparent implementation of the published GEMS distance-weighted Tversky index."""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Any

import numpy as np
from scipy import ndimage


@dataclass(frozen=True)
class MetricComponents:
    tp: float
    fp: float
    fn: float
    n_truth: int
    score: float

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _prediction_array(prediction: np.ndarray, valid: np.ndarray | None = None) -> np.ndarray:
    pred = np.asarray(prediction, dtype=np.float32)
    if pred.ndim != 2:
        raise ValueError(f"prediction must be 2-D, got {pred.shape}")
    if valid is not None:
        valid = np.asarray(valid, dtype=bool)
        if valid.shape != pred.shape:
            raise ValueError("valid mask and prediction must have the same shape")
    # A scorer input is expected to be finite and in [0, 1]. For evaluation of
    # legacy rasters only, NaN/nodata outside the domain is treated as zero.
    pred = np.nan_to_num(pred, nan=0.0, posinf=0.0, neginf=0.0, copy=True)
    pred = np.clip(pred, 0.0, 1.0)
    if valid is not None:
        pred[~valid] = 0.0
    return pred


def _truth_array(truth: np.ndarray, shape: tuple[int, int], valid: np.ndarray | None = None) -> np.ndarray:
    g = np.asarray(truth, dtype=bool)
    if g.ndim != 2 or g.shape != shape:
        raise ValueError(f"truth must have shape {shape}, got {g.shape}")
    if valid is not None:
        valid = np.asarray(valid, dtype=bool)
        if valid.shape != shape:
            raise ValueError("valid mask and truth must have the same shape")
        g = g & valid
    return g


def true_positive_credit(pred: np.ndarray, truth: np.ndarray, radius_px: float = 3.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return per-truth-pixel max credit, and the truth row/column coordinates.

    This follows the published max over predictions within the Euclidean-radius
    neighbourhood; it is not a binary dilation or a gradient-threshold proxy.
    """
    if radius_px <= 0:
        raise ValueError("radius_px must be positive")
    rows, cols = np.nonzero(truth)
    best = np.zeros(rows.size, dtype=np.float32)
    if rows.size == 0:
        return best, rows, cols
    height, width = pred.shape
    rmax = int(math.floor(radius_px))
    for dy in range(-rmax, rmax + 1):
        for dx in range(-rmax, rmax + 1):
            distance = math.hypot(dx, dy)
            if distance > radius_px:
                continue
            rr = rows + dy
            cc = cols + dx
            inside = (rr >= 0) & (rr < height) & (cc >= 0) & (cc < width)
            if not inside.any():
                continue
            contribution = pred[rr[inside], cc[inside]] * (1.0 - distance / radius_px)
            best[inside] = np.maximum(best[inside], contribution)
    return best, rows, cols


def metric_components(
    prediction: np.ndarray,
    truth: np.ndarray,
    *,
    valid: np.ndarray | None = None,
    alpha: float = 0.2,
    beta: float = 0.8,
    radius_px: float = 3.0,
    epsilon: float = 1e-12,
) -> MetricComponents:
    """Compute official DTI components on one raster / scored domain."""
    pred = _prediction_array(prediction, valid)
    g = _truth_array(truth, pred.shape, valid)
    credit, _, _ = true_positive_credit(pred, g, radius_px)
    tp = float(np.sum(credit, dtype=np.float64))
    n_truth = int(np.count_nonzero(g))
    fn = float(n_truth - tp)
    if n_truth:
        distance = ndimage.distance_transform_edt(~g)
        nearest_kernel = np.maximum(1.0 - distance / radius_px, 0.0).astype(np.float32)
    else:
        nearest_kernel = np.zeros(g.shape, dtype=np.float32)
    fp = float(np.sum(pred * (1.0 - nearest_kernel), dtype=np.float64))
    score = tp / (tp + alpha * fp + beta * fn + epsilon) if n_truth else 0.0
    return MetricComponents(tp=tp, fp=fp, fn=fn, n_truth=n_truth, score=float(score))


def spatial_block_components(
    prediction: np.ndarray,
    truth: np.ndarray,
    valid: np.ndarray,
    *,
    n_rows: int = 4,
    n_cols: int = 6,
    guard: int = 3,
    alpha: float = 0.2,
    beta: float = 0.8,
    radius_px: float = 3.0,
    epsilon: float = 1e-12,
) -> tuple[MetricComponents, list[dict[str, Any]]]:
    """Score a frozen spatially blocked holdout and return pooled + per-block DTI.

    A three-cell guard is removed on each *internal* block edge, exactly as the
    project preregistration specifies. Predictions and truth outside the scored
    block cores are ignored. This measures the same published metric on the
    union of disjoint block interiors and avoids cross-block credit transfer.
    """
    pred0 = _prediction_array(prediction, valid)
    g0 = _truth_array(truth, pred0.shape, valid)
    h, w = pred0.shape
    y_edges = np.linspace(0, h, n_rows + 1, dtype=int)
    x_edges = np.linspace(0, w, n_cols + 1, dtype=int)
    scored = np.zeros((h, w), dtype=bool)
    block_masks: list[tuple[int, int, slice, slice]] = []
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
            block_masks.append((iy, ix, rs, cs))

    pred = pred0.copy()
    pred[~scored] = 0.0
    g = g0 & scored
    all_credit, truth_rows, truth_cols = true_positive_credit(pred, g, radius_px)
    n_truth = int(truth_rows.size)
    if n_truth:
        distance = ndimage.distance_transform_edt(~g)
        nearest_kernel = np.maximum(1.0 - distance / radius_px, 0.0).astype(np.float32)
    else:
        nearest_kernel = np.zeros(g.shape, dtype=np.float32)
    fp_map = pred * (1.0 - nearest_kernel)

    tp_by = np.zeros((n_rows, n_cols), dtype=np.float64)
    fn_by = np.zeros((n_rows, n_cols), dtype=np.float64)
    fp_by = np.zeros((n_rows, n_cols), dtype=np.float64)
    if n_truth:
        brow = np.searchsorted(y_edges[1:], truth_rows, side="right")
        bcol = np.searchsorted(x_edges[1:], truth_cols, side="right")
        flatblock = brow * n_cols + bcol
        tp_by.flat[:] = np.bincount(flatblock, weights=all_credit, minlength=n_rows * n_cols)
        counts = np.bincount(flatblock, minlength=n_rows * n_cols)
        fn_by.flat[:] = counts - tp_by.ravel()
    for iy, ix, rs, cs in block_masks:
        fp_by[iy, ix] = np.sum(fp_map[rs, cs], dtype=np.float64)

    blocks: list[dict[str, Any]] = []
    for iy, ix, _, _ in block_masks:
        tp, fp, fn = float(tp_by[iy, ix]), float(fp_by[iy, ix]), float(fn_by[iy, ix])
        den = tp + alpha * fp + beta * fn
        blocks.append({
            "row": iy,
            "col": ix,
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "n_truth": int(round(tp + fn)),
            "score": float(tp / (den + epsilon)) if den > 0 else 0.0,
        })
    total_tp = float(tp_by.sum())
    total_fp = float(fp_by.sum())
    total_fn = float(fn_by.sum())
    den = total_tp + alpha * total_fp + beta * total_fn
    pooled = MetricComponents(
        tp=total_tp,
        fp=total_fp,
        fn=total_fn,
        n_truth=int(round(total_tp + total_fn)),
        score=float(total_tp / (den + epsilon)) if den > 0 else 0.0,
    )
    return pooled, blocks
