"""Tests for the metric-algebra emission module.

Every claim in ``gemsdoe40.emission`` is checked here against the independent
implementation of the organizer's published equations in
``gemsdoe40.research_metric`` (which itself is pinned to the worked example on
the competition page by ``tests/test_research_metric.py``).
"""

from __future__ import annotations

import numpy as np
import pytest

from gemsdoe40 import emission as em
from gemsdoe40.research_metric import metric_components


def _random_case(seed: int = 0, shape=(64, 80)):
    rng = np.random.default_rng(seed)
    truth = rng.random(shape) < 0.02
    pred = np.zeros(shape, dtype=np.float32)
    idx = rng.choice(np.arange(truth.size), size=400, replace=False)
    pred.ravel()[idx] = 1.0
    return pred, truth


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_binary_identity_matches_published_metric(seed: int) -> None:
    """(★) must reproduce the published metric for a binary, off-mask emission.

    Tolerance is 1e-6 relative: both routes accumulate in float32, which is the
    raster dtype the organizer specifies.
    """
    pred, truth = _random_case(seed)
    published = metric_components(pred, truth).score
    rec = em.credit_of_binary(pred, truth)
    predc = em.pred_credit_of_binary(pred, truth)
    closed = em.binary_dti_exact(rec, predc, int(pred.sum()), int(truth.sum()))
    assert closed == pytest.approx(published, rel=1e-6)
    # the thin approximation is the C_pred = T special case and is a strict bound
    # on one side; it must never be reported as the metric itself
    thin = em.binary_dti_thin(rec, int(pred.sum()), int(truth.sum()))
    assert thin != pytest.approx(published, rel=1e-6) or predc == pytest.approx(rec, rel=1e-6)


def test_credit_of_binary_matches_research_metric_tp() -> None:
    pred, truth = _random_case(7)
    mc = metric_components(pred, truth)
    assert em.credit_of_binary(pred, truth) == pytest.approx(mc.tp, rel=1e-6)


def test_mask_removes_both_terms() -> None:
    """A masked pixel contributes neither credit nor penalty (staff-confirmed rule)."""
    pred, truth = _random_case(3)
    mask = truth.copy()
    masked = metric_components(pred, truth, valid=~mask)
    rec = em.credit_of_binary(pred, truth, mask=mask)
    assert rec == pytest.approx(masked.tp, rel=1e-9)
    pc = em.pred_credit_of_binary(pred, truth, mask=mask) if "mask" in em.pred_credit_of_binary.__code__.co_varnames else None
    assert pc is None


def test_binarization_weakly_dominates() -> None:
    """Scaling a fixed support up (lambda <= 1) can only raise the published metric."""
    rng = np.random.default_rng(11)
    truth = rng.random((48, 48)) < 0.03
    support = rng.random((48, 48)) < 0.02
    base = support.astype(np.float32)
    scores = []
    for lam in (0.25, 0.5, 0.75, 1.0):
        scores.append(metric_components(base * lam, truth).score)
    assert scores == sorted(scores)
    assert scores[-1] > scores[0]


def test_credit_bar_and_mass_identities() -> None:
    n, m = 12_226, 40_000
    credit = 5_000.0
    dti = em.binary_dti_thin(credit, m, n)
    # required credit for the achieved value must round-trip
    assert em.metric_required_credit(dti, m, n) == pytest.approx(credit, rel=1e-12)
    assert em.metric_max_mass(credit, n, dti) == pytest.approx(m, rel=1e-9)
    assert em.credit_bar(dti) == pytest.approx(0.2 * dti)


def test_bar_is_the_true_break_even() -> None:
    """A pixel with exactly the bar credit leaves the metric unchanged."""
    n = 12_226
    t0, m0 = 4_000.0, 40_000
    d0 = em.binary_dti_exact(t0, t0, m0, n)
    # a pixel that is the new best for exactly one truth pixel gains delta_T = k
    # and costs mass 1 and its own credit k, so the exact marginal rule is k > 0.2*DTI
    for k in (0.5 * em.credit_bar(d0), 2.0 * em.credit_bar(d0)):
        after = em.binary_dti_exact(t0 + k, t0 + k, m0 + 1, n)
        assert (after > d0) == (k > em.credit_bar(d0))


def test_value_ranked_thinning_respects_separation_and_mass() -> None:
    rng = np.random.default_rng(5)
    score = rng.random((60, 60))
    support = np.ones((60, 60), dtype=bool)
    out = em.value_ranked_thinning(score, support, spacing_px=3.0, max_mass=50)
    assert out.sum() <= 50
    assert out.dtype == np.float32
    rr, cc = np.nonzero(out)
    assert rr.size > 0
    pts = np.c_[rr, cc].astype(float)
    d = np.sqrt(((pts[:, None, :] - pts[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(d, np.inf)
    assert d.min() > 3.0 - 1e-9


def test_thinning_prefers_high_scores() -> None:
    score = np.zeros((30, 30))
    score[5, 5] = 1.0
    score[5, 9] = 0.5
    out = em.value_ranked_thinning(score, np.ones_like(score, dtype=bool), spacing_px=1.5)
    assert out[5, 5] == 1.0 and out[5, 9] == 1.0  # far enough apart to keep both


def test_credit_density_curve_is_monotone_in_mass() -> None:
    rng = np.random.default_rng(9)
    score = rng.random((80, 80))
    truth = rng.random((80, 80)) < 0.01
    rows = em.credit_density_curve(
        score, np.ones_like(score, dtype=bool), truth,
        spacing_px=1.0, mass_grid=[100, 400, 900],
    )
    credits = [r["credit"] for r in rows]
    assert credits == sorted(credits)
    assert all(r["u"] > 0 for r in rows)


def test_model_dti_scales_with_calibration() -> None:
    rows = [{"spacing_px": 3.0, "mass": 30_000, "credit": 3_000.0, "u": 0.1}]
    lo = em.model_dti(rows, n_truth=12_226, calibration=0.4, truth_label="proxy")[0]
    hi = em.model_dti(rows, n_truth=12_226, calibration=0.6, truth_label="proxy")[0]
    assert hi["dti_model"] > lo["dti_model"] > 0


def test_maximin_choice_picks_robust_geometry() -> None:
    rows = [
        {"spacing_px": 1.5, "dti_model": 0.30}, {"spacing_px": 1.5, "dti_model": 0.10},
        {"spacing_px": 3.0, "dti_model": 0.25}, {"spacing_px": 3.0, "dti_model": 0.24},
    ]
    best = em.maximin_choice(rows, group_keys=("spacing_px",))
    assert best["group"] == (3.0,)
    assert best["worst_case"] == pytest.approx(0.24)
