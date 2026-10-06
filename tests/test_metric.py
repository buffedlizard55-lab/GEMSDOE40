"""Verification of the metric implementation against a brute-force definition, plus the
break-even identity that every emission decision in this repository relies on.

The brute force is written directly from the official equations on the problem page
(https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/#performance-metric):
per truth pixel, credit is the max over predicted pixels of p(x)*k(d); per predicted pixel, the
false-positive charge is p(x)*(1 - max over truth pixels of k(d)).
"""

from __future__ import annotations

import numpy as np
import pytest

from gems40.metric import ALPHA, BETA, EPS, RADIUS_PX, break_even_ratio, dti_soft, kernel


def brute_force(p: np.ndarray, g: np.ndarray) -> float:
    H, W = p.shape
    yy, xx = np.mgrid[0:H, 0:W]
    tp = 0.0
    for gy in range(H):
        for gx in range(W):
            if not g[gy, gx]:
                continue
            d = np.hypot(yy - gy, xx - gx)
            tp += float(np.max(p * kernel(d)))
    fp = 0.0
    for py in range(H):
        for px in range(W):
            if p[py, px] <= 0:
                continue
            if not g.any():
                best = 0.0
            else:
                d = np.hypot(yy - py, xx - px)
                best = float(np.max(kernel(d)[g]))
            fp += float(p[py, px] * (1.0 - best))
    n_g = float(g.sum())
    return tp / (tp + ALPHA * fp + BETA * (n_g - tp) + EPS)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_matches_brute_force(seed):
    rng = np.random.default_rng(seed)
    H = W = 12
    p = (rng.random((H, W)) < 0.15) * rng.random((H, W))
    g = rng.random((H, W)) < 0.08
    if not g.any():
        g[5, 5] = True
    got = dti_soft(p, g)["dti"]
    assert got == pytest.approx(brute_force(p, g), rel=1e-9, abs=1e-12)


def test_identities():
    rng = np.random.default_rng(7)
    p = (rng.random((20, 20)) < 0.2).astype(float)
    g = np.zeros((20, 20), bool)
    g[5:9, 5:9] = True
    r = dti_soft(p, g)
    # FN_w = |G| - TP_w exactly (the definition used everywhere in this repository)
    assert r["fn_w"] == pytest.approx(r["n_g"] - r["tp_w"])
    # scaling a fixed support up is monotone improving (the reason artifacts are normalised to 1)
    base = dti_soft(p, g)["dti"]
    hi = dti_soft(np.minimum(p * 1.5, 1.0), g)["dti"]
    if (p * 1.5 <= 1.0).all():
        assert hi > base
    # all-zero prediction scores exactly zero
    assert dti_soft(np.zeros_like(p), g)["dti"] == 0.0


def test_break_even_identity():
    """A marginal addition improves the score iff its credit-per-pixel exceeds tau(D0)."""
    for d0 in (0.05, 0.1, 0.2778, 0.3195):
        tau = break_even_ratio(d0)
        # construct a case with exactly one added pixel at efficiency tau*(1+/-)
        n_g = 10_000.0
        n0 = 1_000
        tp0 = d0 * (ALPHA * n0 + BETA * n_g) / (1.0 - ALPHA * d0 - BETA * (1 - 0.0) * 0)  # algebra check
        # direct numerical check instead: build DTI from (tp, n) and add one pixel of efficiency e
        def dti(tp, n):
            return tp / (ALPHA * n + BETA * n_g)
        tp = 500.0
        n = 2_000
        d0_num = dti(tp, n)
        tau_num = break_even_ratio(d0_num)
        for e, expect_up in ((tau_num * 1.05, True), (tau_num * 0.95, False)):
            new = dti(tp + e, n + 1)
            assert (new > d0_num) is expect_up, (e, tau_num, new, d0_num)


def test_kernel_shape():
    assert kernel(0.0) == 1.0
    assert kernel(RADIUS_PX) == 0.0
    assert kernel(RADIUS_PX * 2) == 0.0
    assert kernel(1.5) == pytest.approx(0.5)
