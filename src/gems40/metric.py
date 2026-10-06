"""Official Distance-Weighted Tversky Index (DTI) for the DOE GEMS Prize (DrivenData #306).

Source of the definition (verified line-by-line against the two official pages):
  * problem description, "Performance metric" section
    https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/
  * organizer clarification that **mapped (knowable) faults are masked out of scoring**
    https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516

Verbatim behaviour implemented here (identical to the organiser's published equations, and
cross-checked against the two independent re-implementations vendored by sibling repos
``5GEMSDOE/scripts/vendor/gems_eval/dti.py`` and ``GEMSDOE32/src/gems32/metric.py``):

    k(d)  = max(1 - d / R, 0)                     R = 300 m = 3 px at the 100 m grid
    TP_w  = sum over truth pixels g of  max over predicted pixels x of  p(x) * k(d(x,g))
    FP_w  = sum over predicted pixels x with p(x) > 0 of  p(x) * (1 - max over g of k(d(x,g)))
    FN_w  = |G| - TP_w
    DTI   = TP_w / (TP_w + alpha*FP_w + beta*FN_w + eps)
          = TP_w / (alpha*(TP_w + FP_w) + beta*|G| + eps)        with alpha = 0.2, beta = 0.8

Notes that drive every emission decision in this repository:
  1. ``TP_w`` takes a **max** over predictions for each truth pixel, while ``FP_w`` **sums**
     over predictions.  Redundant mass near truth is therefore a discount, not a penalty.
  2. ``alpha``/``beta`` are 0.2/0.8, i.e. a missing truth pixel costs 4x an equal false-positive
     pixel -> recall-leaning *if* the mass lands near real structure.
  3. Because every term is linear in ``p`` for a fixed support and ``beta*|G|`` is a constant,
     for a fixed spatial pattern scaling all values up is strictly monotone improving
     (d DTI / d lambda > 0).  The metric optimum is therefore at the top of the [0, 1] range,
     which is why the shipped artifact is normalised with max = 1.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import distance_transform_edt

ALPHA = 0.2
BETA = 0.8
RADIUS_PX = 3.0
EPS = 1e-12

# 5x5x... offsets with k(d) > 0 for R = 3 px (25 px disk-ish stencil, exactly the official one)
_OFFSETS: tuple[tuple[int, int, float], ...] = tuple(
    (dy, dx, max(1.0 - float(np.hypot(dy, dx)) / RADIUS_PX, 0.0))
    for dy in range(-3, 4)
    for dx in range(-3, 4)
    if max(1.0 - float(np.hypot(dy, dx)) / RADIUS_PX, 0.0) > 0.0
)


def kernel(d):
    """Triangular kernel k(d) = max(1 - d / 3 px, 0)."""
    return np.maximum(1.0 - np.asarray(d, dtype=np.float64) / RADIUS_PX, 0.0)


def dti_soft(pred: np.ndarray, truth: np.ndarray, valid: np.ndarray | None = None,
             known: np.ndarray | None = None, alpha: float = ALPHA, beta: float = BETA) -> dict:
    """Exact official DTI for a *soft* prediction raster (values in [0, 1])."""
    pred = np.asarray(pred, dtype=np.float64)
    truth = np.asarray(truth)
    valid = np.ones(pred.shape, bool) if valid is None else np.asarray(valid, bool)
    known = np.zeros(pred.shape, bool) if known is None else np.asarray(known, bool)
    active = valid & ~known
    p = np.where(active & np.isfinite(pred), pred, 0.0)
    g = active & (np.asarray(truth) > 0)
    n_g = int(g.sum())
    if n_g == 0:
        return dict(dti=0.0, tp_w=0.0, fp_w=float(p.sum()), n_g=0, n_pred_mass=float(p.sum()))
    yy, xx = np.nonzero(g)
    credit = np.zeros(yy.size, dtype=np.float64)
    H, W = p.shape
    for dy, dx, k in _OFFSETS:
        ny, nx = yy + dy, xx + dx
        ok = (ny >= 0) & (ny < H) & (nx >= 0) & (nx < W)
        credit[ok] = np.maximum(credit[ok], p[ny[ok], nx[ok]] * k)
    tp_w = float(credit.sum())
    # FP mass = every prediction pixel weighted by (1 - best kernel value to any truth pixel)
    d_g = distance_transform_edt(~g)
    near = kernel(d_g)                     # per prediction pixel: max_g k(d)
    fp_w = float((p * (1.0 - near)).sum())
    fn_w = float(n_g) - tp_w
    dti = tp_w / (tp_w + alpha * fp_w + beta * fn_w + EPS)
    return dict(dti=float(dti), tp_w=tp_w, fp_w=fp_w, fn_w=fn_w, n_g=n_g,
                n_pred_mass=float(p.sum()), n_pred_px=int((p > 0).sum()))


def dti_binary(pred_mask: np.ndarray, truth: np.ndarray) -> dict:
    """Official DTI for a binary support mask (equivalent to dti_soft with p in {0,1})."""
    return dti_soft(np.asarray(pred_mask, dtype=np.float64), truth)


def break_even_ratio(dti0: float, alpha: float = ALPHA) -> float:
    """tau = alpha*D0 / (1 - alpha*D0): the credit/FP ratio above which adding mass helps.

    Derived from d DTI > 0 for a marginal addition (Delta TP_w, Delta FP_w) at operating point D0.
    """
    return float(alpha * dti0 / (1.0 - alpha * dti0))
