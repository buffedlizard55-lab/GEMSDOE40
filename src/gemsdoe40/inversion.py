#!/usr/bin/env python3
"""Invert the published leaderboard scores into a spatial truth-density estimate.

Why this exists
---------------
Each scored prediction's pixels are byte-recoverable from the pinned public corpus and
the organizer metric is a known function of the hidden truth:

    DTI = T / (0.2*T + 0.2*(M - C) + 0.8*N)
    T = credit arriving at truth pixels from the prediction (<= M)
    C = credit seen by the emitted dots from the nearest truth pixel
    M = emitted mass, N = number of truth pixels in the scored domain

so a set of scored predictions are *measurements of the hidden truth*.  Modelling the
hidden truth as an independent-pixel probability field ``rho`` on a coarse block grid:

    T_k = sum_b rho_b * A_kb                    (A = block sums of the kernel field)
    C_k = sum_b (1 - exp(-rho_b * s)) * D_kb    (D = anchor dot counts per block)
    N   = sum_b rho_b * m_b                     (m = scoreable pixels per block)
    s   = kernel mass of the 3 px disk (~9.42)

Both ``T`` and ``C`` account for the same pixels, so ``M - C`` is the expected
false-positive mass.  The map is solved with non-negativity, Laplacian smoothing and a
weak prior toward an independent public fault map, then validated by leave-one-out
prediction of the published scores.

Nothing here is a competition score: the model is an instrument, and it is only usable
if its leave-one-out ranking of the scored anchors is good.
"""
from __future__ import annotations

import numpy as np

RADIUS_PX = 3.0
ALPHA = 0.2
BETA = 0.8
RHO_MAX = 0.5


def kernel_mass(radius_px: float = RADIUS_PX) -> float:
    """Sum of the metric's triangular kernel over the integer offsets it supports."""
    rad = int(np.floor(radius_px))
    total = 0.0
    for dy in range(-rad, rad + 1):
        for dx in range(-rad, rad + 1):
            d = float(np.hypot(dx, dy))
            if d <= radius_px:
                total += 1.0 - d / radius_px
    return total


class InversionModel:
    """Block-grid truth-density model of the organizer metric."""

    def __init__(self, block: int, shape: tuple[int, int]):
        self.block = int(block)
        self.shape = tuple(int(v) for v in shape)
        self.rows = int(np.ceil(self.shape[0] / self.block))
        self.cols = int(np.ceil(self.shape[1] / self.block))
        self.s = kernel_mass()
        self.n_blocks = self.rows * self.cols
        self._idx = None
        self.anchors: list[dict] = []
        self.m_b = np.zeros(self.n_blocks)
        self.prior_density = None

    def block_index(self) -> np.ndarray:
        if self._idx is None:
            rr, cc = np.indices(self.shape)
            self._idx = ((rr // self.block) * self.cols + (cc // self.block)).astype(np.int64)
        return self._idx

    def add_anchor(self, name: str, dots: np.ndarray, kernel_field: np.ndarray,
                   live_score: float) -> None:
        flat = self.block_index().ravel()
        self.anchors.append({
            "name": name,
            "live_score": float(live_score),
            "mass": int(dots.sum()),
            "dot_counts": np.bincount(flat[dots.ravel()], minlength=self.n_blocks).astype(np.float64),
            "kernel_sum": np.bincount(flat, weights=kernel_field.ravel(),
                                      minlength=self.n_blocks).astype(np.float64),
        })

    def set_scoreable(self, scoreable: np.ndarray) -> None:
        flat = self.block_index().ravel()
        self.m_b = np.bincount(flat[scoreable.ravel()], minlength=self.n_blocks).astype(np.float64)

    def set_prior(self, prior_density: np.ndarray) -> None:
        """Per-block prior density on the block grid (see :meth:`density_of_mask`)."""
        arr = np.asarray(prior_density, dtype=np.float64)
        if arr.shape != (self.rows, self.cols):
            raise ValueError(f"prior density must be {(self.rows, self.cols)}, got {arr.shape}")
        self.prior_density = arr

    def density_of_mask(self, mask: np.ndarray) -> np.ndarray:
        """Fraction of scoreable pixels inside each block covered by ``mask``."""
        flat = self.block_index().ravel()
        counts = np.bincount(flat[np.asarray(mask, bool).ravel()],
                             minlength=self.n_blocks).astype(np.float64)
        return (counts / np.maximum(self.m_b, 1.0)).reshape(self.rows, self.cols)

    # ---- forward model --------------------------------------------------- #
    def score(self, rho: np.ndarray, anchor: dict) -> float:
        rho = np.asarray(rho, dtype=np.float64)
        if rho.shape != (self.rows, self.cols):
            raise ValueError(f"rho must be {(self.rows, self.cols)}, got {rho.shape}")
        r = rho.ravel()
        t = float((r * anchor["kernel_sum"]).sum())
        c = float(((1.0 - np.exp(-r * self.s)) * anchor["dot_counts"]).sum())
        n = float((r * self.m_b).sum())
        den = ALPHA * t + ALPHA * (float(anchor["mass"]) - c) + BETA * n
        if den <= 1e-9:
            return 0.0
        return min(t / den, 1.0)

    def score_all(self, rho: np.ndarray) -> dict[str, float]:
        return {a["name"]: self.score(rho, a) for a in self.anchors}


def _objective(theta: np.ndarray, model: InversionModel, anchors: list[dict],
               peaks: np.ndarray, lam: float, mu: float, prior: np.ndarray | None):
    shape = (model.rows, model.cols)
    rho = (RHO_MAX / (1.0 + np.exp(-theta))).reshape(shape)
    r = rho.ravel()
    use_prior = prior is not None and mu > 0
    res = np.zeros(len(anchors) + model.n_blocks + (model.n_blocks if use_prior else 0))
    grad = np.zeros(model.n_blocks)

    i = 0
    for anchor, y in zip(anchors, peaks):
        ks, dc = anchor["kernel_sum"], anchor["dot_counts"]
        t = float((r * ks).sum())
        g = 1.0 - np.exp(-r * model.s)
        c = float((g * dc).sum())
        n = float((r * model.m_b).sum())
        m = float(anchor["mass"])
        den = ALPHA * t + ALPHA * (m - c) + BETA * n
        raw = 0.0 if den <= 1e-9 else t / den
        d = min(raw, 1.0)
        res[i] = d - y
        i += 1
        if den > 1e-9 and raw <= 1.0:
            dg = model.s * np.exp(-r * model.s)
            dden = ALPHA * ks - ALPHA * dg * dc + BETA * model.m_b
            grad += (2.0 * (d - y)) * ((ks * den - t * dden) / (den * den))

    # Laplacian smoothing with Neumann boundaries (self-adjoint => grad = lam * lap2).
    padded = np.pad(rho, 1, mode="edge")
    lap = (padded[1:-1, :-2] + padded[1:-1, 2:] + padded[:-2, 1:-1] + padded[2:, 1:-1]
           - 4.0 * rho)
    res[i:i + model.n_blocks] = (np.sqrt(lam) * lap).ravel()
    i += model.n_blocks
    lp = np.pad(lap, 1, mode="edge")
    lap2 = (lp[1:-1, :-2] + lp[1:-1, 2:] + lp[:-2, 1:-1] + lp[2:, 1:-1] - 4.0 * lap)
    grad += lam * lap2.ravel()

    if use_prior:
        diff = rho - prior
        res[i:i + model.n_blocks] = (np.sqrt(mu) * diff).ravel()
        grad += mu * diff.ravel()

    drho = rho * (1.0 - rho / RHO_MAX)
    grad *= drho.ravel()
    return float((res ** 2).sum()), grad


def fit(model: InversionModel, subset: list[str] | None = None, *, lam: float = 1e3,
        mu: float = 0.0, init: np.ndarray | None = None, max_iter: int = 300):
    """Fit the block density by projected L-BFGS on a logit parameterisation."""
    from scipy import optimize

    anchors = [a for a in model.anchors if subset is None or a["name"] in subset]
    peaks = np.array([a["live_score"] for a in anchors], float)
    if init is None:
        init = (model.prior_density.copy() if model.prior_density is not None
                else np.full((model.rows, model.cols), 0.02))
    init = np.clip(init, 1e-6, RHO_MAX - 1e-6)
    theta0 = np.log(init / (RHO_MAX - init)).ravel()
    prior = model.prior_density
    out = optimize.minimize(_objective, theta0, jac=True, method="L-BFGS-B",
                            args=(model, anchors, peaks, lam, mu, prior),
                            options=dict(maxiter=max_iter, maxfun=max_iter * 2))
    rho = (RHO_MAX / (1.0 + np.exp(-out.x))).reshape(model.rows, model.cols)
    return rho, out


def spearman(a, b) -> float | None:
    from scipy import stats

    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) < 3 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return None
    return float(stats.spearmanr(a, b).statistic)


def leave_one_out(model: InversionModel, *, lam: float, mu: float, max_iter: int = 300,
                  init_from_full: bool = True) -> dict:
    names = [a["name"] for a in model.anchors]
    base = None
    if init_from_full:
        base, _ = fit(model, None, lam=lam, mu=mu, max_iter=max_iter)
    preds, actual = [], []
    for held in names:
        rho, _ = fit(model, [n for n in names if n != held], lam=lam, mu=mu,
                     init=base, max_iter=max_iter)
        held_anchor = next(a for a in model.anchors if a["name"] == held)
        preds.append(model.score(rho, held_anchor))
        actual.append(held_anchor["live_score"])
    return {"loo_spearman": spearman(preds, actual),
            "loo_rmse": float(np.sqrt(np.mean((np.array(preds) - np.array(actual)) ** 2))),
            "pred": preds, "actual": actual, "names": names}
