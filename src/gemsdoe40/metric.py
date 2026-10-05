"""Official distance-weighted Tversky index (DTI).

Source (verbatim equations): DrivenData problem description, Page 967
https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/#performance-metric

    k(d) = max(1 - d/R, 0)   R = 300 m = 3 px at 100 m
    TP_w = sum_{g in G} max_{x: d(x,g)<=R} p(x) k(d(x,g))
    FP_w = sum_{x: p(x)>0} p(x) [1 - max_g k(d(x,g))]
    FN_w = |G| - TP_w
    DTI(α=0.2, β=0.8) = TP_w / (TP_w + α FP_w + β FN_w)

Worked example on that page: TP=3.00, FP=1.89, FN=2.00 → DTI = 0.60
(reproduced in tests/test_metric.py).
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import distance_transform_edt

ALPHA = 0.2
BETA = 0.8
RADIUS_PX = 3.0
EPS = 1e-12


def kernel(d, radius: float = RADIUS_PX):
    return np.maximum(1.0 - np.asarray(d, dtype=np.float64) / radius, 0.0)


def dti_from_components(tp: float, fp: float, fn: float,
                        alpha: float = ALPHA, beta: float = BETA) -> float:
    return float(tp / (tp + alpha * fp + beta * fn + EPS))


def dti_exact(pred, truth, valid=None, alpha: float = ALPHA, beta: float = BETA) -> dict:
    """Exact DTI for soft or binary predictions in [0, 1]."""
    p = np.asarray(pred, dtype=np.float64)
    g = np.asarray(truth, dtype=bool)
    if p.shape != g.shape:
        raise ValueError("prediction and truth shape mismatch")
    if valid is not None:
        valid = np.asarray(valid, dtype=bool)
        p = np.where(valid, p, 0.0)
        g = g & valid
    if np.isnan(p).any() or (p < 0).any() or (p > 1).any():
        # NaNs outside the scored domain should already have been zeroed.
        bad = ~np.isfinite(p) | (p < 0) | (p > 1)
        if bad.any():
            raise ValueError("predictions must be finite and in [0, 1] in the scored domain")
    n = int(g.sum())
    if n == 0:
        return dict(tp=0.0, fp=float(p.sum()), fn=0.0, n_truth=0, dti=0.0)
    r = int(np.ceil(RADIUS_PX))
    yy, xx = np.nonzero(g)
    credit = np.zeros(n, dtype=np.float64)
    H, W = p.shape
    for dy in range(-r, r + 1):
        for dx in range(-r, r + 1):
            k = float(kernel(np.hypot(dy, dx)))
            if k <= 0:
                continue
            ny, nx = yy + dy, xx + dx
            ok = (ny >= 0) & (ny < H) & (nx >= 0) & (nx < W)
            credit[ok] = np.maximum(credit[ok], p[ny[ok], nx[ok]] * k)
    tp = float(credit.sum())
    fn = float(n) - tp
    dist = distance_transform_edt(~g)
    fp = float((p * (1.0 - kernel(dist))).sum())
    return dict(tp=tp, fp=fp, fn=fn, n_truth=n, dti=dti_from_components(tp, fp, fn, alpha, beta))


def dti_binary(pred_bool, truth, valid=None) -> dict:
    """Fast exact DTI when predictions are {0, 1}."""
    p = np.asarray(pred_bool, dtype=bool)
    g = np.asarray(truth, dtype=bool)
    if valid is not None:
        valid = np.asarray(valid, dtype=bool)
        p = p & valid
        g = g & valid
    n = int(g.sum())
    if n == 0:
        return dict(tp=0.0, fp=float(p.sum()), fn=0.0, n_truth=0, dti=0.0)
    if not p.any():
        return dict(tp=0.0, fp=0.0, fn=float(n), n_truth=n, dti=0.0)
    tp = float(kernel(distance_transform_edt(~p)[g]).sum())
    fn = float(n) - tp
    fp = float((1.0 - kernel(distance_transform_edt(~g)[p])).sum())
    return dict(tp=tp, fp=fp, fn=fn, n_truth=n, dti=dti_from_components(tp, fp, fn))


def marginal_bar(current_dti: float, alpha: float = ALPHA) -> float:
    """Minimum kernel credit k for an added unit of mass to raise DTI.

    dDTI > 0  <=>  k > α · DTI   (denominator rises by exactly α regardless of k).
    """
    return alpha * current_dti
