"""Calibrate a surrogate-credit measurement to the official leaderboard scores.

The organiser's metric is a ratio in which the two unknown quantities are the
hidden label mass ``G`` (pixels of expert-added faults inside the scored
chunk) and how much kernel credit a candidate's pixels actually earn on that
hidden set.  A surrogate truth set cannot tell us either, but it *can* rank
placements, because the surrogate credit is a monotone function of the hidden
credit for structures of the same geometry.

So the model has exactly two free parameters:

    TP_hidden(N)  = gamma * N * h_surrogate(N)
    FP_hidden(N)  = N * (1 - gamma * kappa_surrogate(N))
    DTI(N)        = 5 TP / (TP + FP + 4 G)

where ``h_surrogate`` is the mean true-positive credit per emitted dot and
``kappa_surrogate`` the mean own-credit per emitted dot (both measured on the
surrogate truth), ``gamma`` maps surrogate credit to hidden credit, and ``G``
is the hidden truth pixel count.  Both are fitted by least squares to the
(emitted pixel count, official score) pairs of prior submissions whose scores
are known, excluding candidates built from the surrogate catalogue itself
(circular) and byte-identical duplicates.

Leave-one-out fits are reported so the reader can see how much of the fit is
signal: the model is a calibration, not a promise.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

ALPHA = 0.2
BETA = 0.8


def dti_from_components(tp: float, fp: float, n_truth: float,
                        alpha: float = ALPHA, beta: float = BETA) -> float:
    fn = max(float(n_truth) - float(tp), 0.0)
    denom = tp + alpha * fp + beta * fn
    return float(tp / denom) if denom > 0 else 0.0


def predict_score(n_px: float, h_sur: float, kappa_sur: float, gamma: float, g_hidden: float,
                  alpha: float = ALPHA, beta: float = BETA) -> float:
    """Model DTI for a candidate with ``n_px`` emitted dots."""
    tp = gamma * float(n_px) * float(h_sur)
    fp = float(n_px) * (1.0 - gamma * float(kappa_sur))
    tp = float(np.clip(tp, 0.0, g_hidden))
    fp = max(fp, 0.0)
    return dti_from_components(tp, fp, g_hidden, alpha, beta)


@dataclass
class Fit:
    gamma: float
    g_hidden: float
    rmse: float
    r2: float
    n: int
    loo_rmse: float
    loo_spearman: float
    residuals: list


def _pack(params, obs):
    gamma, g_hidden = params
    if gamma <= 0 or g_hidden <= 0:
        return None
    return np.array([predict_score(o["n_px"], o["h_sur"], o["kappa_sur"], gamma, g_hidden)
                     for o in obs])


def fit_model(obs: list[dict], *, gamma_bounds=(0.05, 3.0), g_bounds=(2e3, 4e5),
              n_grid: int = 40, with_loo: bool = True) -> Fit:
    """Grid + local refinement fit of (gamma, G); returns LOO diagnostics.

    ``with_loo`` is disabled for the internal leave-one-out refits; without that
    guard the recursion is factorial.
    """
    obs = [o for o in obs if o.get("use_for_fit", True)]
    if len(obs) < 3:
        raise ValueError("need at least 3 observations to fit")
    y = np.array([o["official"] for o in obs], dtype=float)
    gammas = np.geomspace(gamma_bounds[0], gamma_bounds[1], n_grid)
    gs = np.geomspace(g_bounds[0], g_bounds[1], n_grid * 3)

    def loss(gamma, g):
        pred = np.array([predict_score(o["n_px"], o["h_sur"], o["kappa_sur"], gamma, g) for o in obs])
        return float(np.mean((pred - y) ** 2))

    best = (np.inf, gammas[0], gs[0])
    for gamma in gammas:
        for g in gs:
            l = loss(gamma, g)
            if l < best[0]:
                best = (l, gamma, g)
    # local refinement
    _, gamma0, g0 = best
    for _ in range(4):
        gammas = np.linspace(gamma0 * 0.7, gamma0 * 1.4, 15)
        gs = np.linspace(g0 * 0.7, g0 * 1.4, 15)
        for gamma in gammas:
            for g in gs:
                l = loss(gamma, g)
                if l < best[0]:
                    best = (l, gamma, g)
        _, gamma0, g0 = best
    rmse, gamma, g_hidden = best
    pred = np.array([predict_score(o["n_px"], o["h_sur"], o["kappa_sur"], gamma, g_hidden) for o in obs])
    ss_res = float(np.sum((pred - y) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    # leave-one-out
    loo_pred, loo_true = [], []
    if with_loo and len(obs) >= 4:
        for i in range(len(obs)):
            sub = obs[:i] + obs[i + 1:]
            try:
                f = fit_model(sub, gamma_bounds=gamma_bounds, g_bounds=g_bounds, n_grid=20,
                              with_loo=False)
            except ValueError:
                continue
            loo_pred.append(predict_score(obs[i]["n_px"], obs[i]["h_sur"], obs[i]["kappa_sur"],
                                          f.gamma, f.g_hidden))
            loo_true.append(obs[i]["official"])
    loo_rmse = float(np.sqrt(np.mean((np.array(loo_pred) - np.array(loo_true)) ** 2))) if loo_pred else float("nan")
    loo_spearman = float("nan")
    if len(loo_pred) >= 3:
        from scipy.stats import spearmanr
        loo_spearman = float(spearmanr(loo_pred, loo_true).statistic)
    return Fit(gamma=float(gamma), g_hidden=float(g_hidden), rmse=float(np.sqrt(rmse)),
               r2=float(r2), n=len(obs), loo_rmse=loo_rmse, loo_spearman=loo_spearman,
               residuals=[dict(id=o.get("id"), official=o["official"], modelled=float(p),
                               n_px=o["n_px"], h_sur=o["h_sur"], kappa_sur=o["kappa_sur"])
                          for o, p in zip(obs, pred)])


def optimal_budget(curve: list[dict], gamma: float, g_hidden: float) -> dict:
    """Given (n_px, h_sur, kappa_sur) points along the ranked dot list, best modelled score."""
    best = None
    for c in curve:
        s = predict_score(c["n_px"], c["h_sur"], c["kappa_sur"], gamma, g_hidden)
        if best is None or s > best["modelled_dti"]:
            best = {**{k: v for k, v in c.items() if k != "modelled_dti"}, "modelled_dti": s}
    margin = None
    if best is not None:  # margin vs the next-best distinct budget of the same curve
        others = sorted((predict_score(c["n_px"], c["h_sur"], c["kappa_sur"], gamma, g_hidden)
                         for c in curve), reverse=True)
        margin = float(others[0] - others[1]) if len(others) > 1 else None
    if best is not None:
        best["margin_over_next"] = margin
    return best or {}
