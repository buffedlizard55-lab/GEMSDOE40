"""Spatially blocked logistic discriminant, used as the H41 emission ranker.

Why a discriminant is allowed here
----------------------------------
The competition rules permit any additional data source whose licence allows use in the
challenge ("External datasets", problem description page 967).  The USGS State Geologic
Map Compilation (SGMC) fault layer is public-domain USGS material, it is *not* the
organizer's label raster, and the task's own target is faults that are absent from the
provided USGS catalogue.  Faults in SGMC that are missing from the provided catalogue
are therefore a legitimate, independent, *public* stand-in for the unmapped-fault class.

Two honesty requirements are built into this module:

* **Spatial blocking.**  Folds are contiguous geographic blocks, never random pixels, so
  a pixel's own fault line and its immediate neighbours are always in the same fold.
  The emitted probability field is out-of-fold everywhere: every block is scored by a
  model that never saw that block.
* **It is not the hidden truth.**  Agreement with SGMC-off-catalogue faults is a
  measurable proxy objective, not the leaderboard.  Every number this module produces is
  labelled as a proxy.

The classifier is deliberately simple (L2-penalised logistic regression fitted by
L-BFGS on a balanced sample), because the audit trail must be reproducible without a
machine-learning framework.
"""
from __future__ import annotations

import numpy as np


def blocked_folds(shape: tuple[int, int], rows: int = 5, cols: int = 6,
                  block_px: int | None = None) -> tuple[np.ndarray, int, int]:
    """Contiguous block fold ids covering ``shape``.

    Returns ``(fold_id, n_rows, n_cols)`` where ``fold_id`` is ``-1`` outside the
    rectangle that the balanced grid covers (the remainder strip is assigned to the
    adjacent block so no pixel is left unscored).
    """
    h, w = shape
    n = rows * cols
    rr, cc = np.indices(shape)
    br = np.minimum((rr * rows) // h, rows - 1)
    bc = np.minimum((cc * cols) // w, cols - 1)
    return (br * cols + bc).astype(np.int32), rows, cols


def fit_logistic(x: np.ndarray, y: np.ndarray, *, l2: float = 1.0, max_iter: int = 500,
                 standardize: bool = True, clip_z: float = 10.0,
                 sample_weight: np.ndarray | None = None) -> dict:
    """L2-penalised logistic regression by L-BFGS with analytic gradients.

    ``x`` is ``(n_samples, n_features)`` float64, ``y`` is 0/1.  Returns the fitted
    coefficients plus the standardisation used, so the same transform can be applied to
    the full grid.
    """
    from scipy import optimize

    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.ndim != 2 or x.shape[0] != y.size:
        raise ValueError("x must be (n, d) and y (n,)")
    # Winsorise heavy-tailed geophysical bands before standardising: a single 1e38
    # sentinel or an extreme spike otherwise dominates the log-loss and overflows it.
    lo = np.percentile(x, 0.5, axis=0)
    hi = np.percentile(x, 99.5, axis=0)
    x = np.clip(x, lo, hi)
    mu = x.mean(axis=0) if standardize else np.zeros(x.shape[1])
    sd = x.std(axis=0) if standardize else np.ones(x.shape[1])
    sd = np.where(sd > 0, sd, 1.0)
    z = np.clip((x - mu) / sd, -clip_z, clip_z)
    z = np.hstack([np.ones((z.shape[0], 1)), z])
    # Class weights keep the balanced sample from being dominated by the majority class;
    # an optional per-sample weight lets the caller pass soft (kernel-credit) labels.
    p1 = max(float((y * (1.0 if sample_weight is None else sample_weight)).sum()
                   / (1.0 if sample_weight is None else sample_weight.sum())), 1e-6)
    wgt = np.where(y > 0.5, 0.5 / p1, 0.5 / max(1.0 - p1, 1e-6))
    if sample_weight is not None:
        wgt = wgt * np.asarray(sample_weight, dtype=np.float64)

    def objective(beta: np.ndarray):
        s = np.clip(z @ beta, -60.0, 60.0)
        # numerically stable log-loss
        ll = np.sum(wgt * (np.logaddexp(0.0, s) - y * s))
        pen = 0.5 * l2 * float(beta[1:] @ beta[1:])
        prob = 1.0 / (1.0 + np.exp(-s))
        grad = z.T @ (wgt * (prob - y))
        grad[1:] += l2 * beta[1:]
        return ll / z.shape[0] + pen / z.shape[0], grad / z.shape[0]

    beta0 = np.zeros(z.shape[1])
    out = optimize.minimize(objective, beta0, jac=True, method="L-BFGS-B",
                            options=dict(maxiter=max_iter))
    return {"beta": out.x, "mu": mu, "sd": sd, "lo": lo, "hi": hi, "clip_z": clip_z,
            "soft_labels": bool(sample_weight is not None),
            "success": bool(out.success),
            "n_iter": int(out.nit), "n_samples": int(z.shape[0]),
            "n_features": int(x.shape[1])}


def predict_logistic(model: dict, x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    lo = model.get("lo")
    hi = model.get("hi")
    if lo is not None and hi is not None:
        x = np.clip(x, lo, hi)
    z = np.clip((x - model["mu"]) / model["sd"], -(model.get("clip_z") or 10.0),
                model.get("clip_z") or 10.0)
    z = np.hstack([np.ones((z.shape[0], 1)), z])
    s = np.clip(z @ model["beta"], -60.0, 60.0)
    return 1.0 / (1.0 + np.exp(-s))


def auc(y: np.ndarray, score: np.ndarray) -> float | None:
    """Rank-based AUC (Mann-Whitney), ties handled by average ranks."""
    from scipy import stats

    y = np.asarray(y, bool)
    if y.all() or not y.any():
        return None
    r = stats.rankdata(score)
    n1 = float(y.sum())
    n0 = float((~y).sum())
    return float((r[y].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def precision_at_top(y: np.ndarray, score: np.ndarray, k: int) -> float:
    y = np.asarray(y, bool)
    k = int(min(k, y.size))
    if k <= 0:
        return 0.0
    idx = np.argsort(-np.asarray(score, float), kind="stable")[:k]
    return float(y[idx].mean())
