"""Spatially blocked holdout for an unsupervised Euler field.

Euler deconvolution does not train on labels, so there is no leakage from
fitting.  The holdout answers a different question: *does the Euler cloud
preferentially fall on held-out catalogue traces more than a same-mass
gradient-threshold field and a same-mass random field?*  A pass licenses
packaging a candidate; it is not a live score.  The hidden test set is
off-catalogue by construction, so a high on-catalogue holdout DTI is only
evidence of *fault-finding skill*, which we then deploy off-catalogue.
"""
from __future__ import annotations

import numpy as np

from .metric import dti_exact, dti_binary


def quadrants(shape: tuple[int, int]) -> dict[str, np.ndarray]:
    h, w = shape
    r = np.arange(h)[:, None]
    c = np.arange(w)[None, :]
    mid_r, mid_c = h // 2, w // 2
    return {
        "NW": (r < mid_r) & (c < mid_c),
        "NE": (r < mid_r) & (c >= mid_c),
        "SW": (r >= mid_r) & (c < mid_c),
        "SE": (r >= mid_r) & (c >= mid_c),
    }


def blocked_dti(pred: np.ndarray, truth: np.ndarray, valid: np.ndarray) -> dict:
    """Four-fold spatial block: score pred against truth inside each quadrant."""
    folds = {}
    dtis = []
    for name, q in quadrants(pred.shape).items():
        mask = valid & q
        res = dti_exact(pred, truth, valid=mask)
        folds[name] = res
        dtis.append(res["dti"])
    return {
        "per_fold": folds,
        "mean_dti": float(np.mean(dtis)),
        "folds_positive_vs_zero": int(sum(d > 0 for d in dtis)),
    }


def same_mass_random(pred: np.ndarray, valid: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Binary random field with the same number of positive pixels as ``pred``."""
    n = int((pred > 0).sum())
    idx = np.flatnonzero(valid.ravel())
    n = min(n, idx.size)
    pick = rng.choice(idx, size=n, replace=False)
    out = np.zeros(pred.shape, dtype=np.float32)
    out.ravel()[pick] = 1.0
    return out


def gradient_baseline(field: np.ndarray, valid: np.ndarray, n_emit: int) -> np.ndarray:
    """Same-count baseline: top-n pixels of |∇field| (the thing Euler is *not*)."""
    f = np.where(valid, field, np.nan)
    gy, gx = np.gradient(np.nan_to_num(f, nan=0.0))
    mag = np.hypot(gx, gy)
    mag = np.where(valid, mag, -np.inf)
    if n_emit <= 0:
        return np.zeros(field.shape, dtype=np.float32)
    flat = mag.ravel()
    # argpartition for the top-n
    n_emit = min(n_emit, int(valid.sum()))
    thresh_idx = np.argpartition(flat, -n_emit)[-n_emit:]
    out = np.zeros(field.shape, dtype=np.float32)
    out.ravel()[thresh_idx] = 1.0
    out = np.where(valid, out, 0.0)
    return out


def evaluate(pred: np.ndarray, truth: np.ndarray, valid: np.ndarray,
             grad_field: np.ndarray | None, seed: int = 40) -> dict:
    rng = np.random.default_rng(seed)
    n_emit = int((pred > 0).sum())
    blocked = blocked_dti(pred, truth, valid)
    rnd = same_mass_random(pred, valid, rng)
    rnd_res = dti_binary(rnd > 0, truth, valid=valid)
    out = {
        "emitted_px": n_emit,
        "blocked": blocked,
        "full_dti": dti_exact(pred, truth, valid=valid),
        "random_same_mass_dti": rnd_res,
        "delta_vs_random": float(blocked["mean_dti"] - rnd_res["dti"]),
    }
    if grad_field is not None and n_emit > 0:
        g = gradient_baseline(grad_field, valid, n_emit)
        g_res = dti_binary(g > 0, truth, valid=valid)
        out["gradient_same_mass_dti"] = g_res
        out["delta_vs_gradient"] = float(out["full_dti"]["dti"] - g_res["dti"])
    return out
