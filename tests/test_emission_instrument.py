"""Regression tests for the live-score-calibrated instrument.

The instrument is fitted outside this repository (sibling repository
``GEMSDOE39``, session 2) to 12 organizer-scored artifacts.  These tests pin
(a) the published anchor predictions, (b) the definition of ``w``, and
(c) the marginal break-even rule, so that any future edit that silently changes
the decision rule fails here rather than in a submission decision.
"""
from __future__ import annotations

import numpy as np
import pytest

from gemsdoe40.emission import (
    calibrated_live_score,
    dot_credit_w,
    instrument_break_even,
)

#: (n_dots, w, published prediction) — values published by the sibling fit.
ANCHORS = [
    (204504, 0.0505, 0.1090),   # blind spacing-5 lattice, live 0.0904
    (166519, 0.0834, 0.1395),   # Hedge-v2 / ens12, live 0.1563
    (161366, 0.0843, 0.1428),   # H25-ctx-ridge, live 0.1280
    (65236, 0.0710, 0.1914),    # H28 dotted ridge, live 0.1839
    (123939, 0.1051, 0.1763),   # H16-1, live 0.1855
    (121131, 0.1035, 0.1782),   # H19-5, live 0.1922
    (60069, 0.0997, 0.2327),    # dotted H19-5 d1.5, live 0.2477
    (44090, 0.0989, 0.2399),    # dotted H19-5 d2.8, live 0.2600
]


@pytest.mark.parametrize("n,w,published", ANCHORS)
def test_published_anchor_predictions_reproduce(n: int, w: float, published: float) -> None:
    assert calibrated_live_score(n, w) == pytest.approx(published, abs=2e-4)


def test_w_is_one_on_truth_and_zero_beyond_the_kernel_radius() -> None:
    truth = np.zeros((21, 21), dtype=bool)
    truth[10, 10] = True
    dots = np.zeros_like(truth)
    dots[10, 10] = True
    assert dot_credit_w(dots, truth) == pytest.approx(1.0)
    far = np.zeros_like(truth)
    far[10, 14] = True                      # 4 px > R = 3 px
    assert dot_credit_w(far, truth) == pytest.approx(0.0)
    mid = np.zeros_like(truth)
    mid[10, 12] = True                      # 2 px -> 1 - 2/3
    assert dot_credit_w(mid, truth) == pytest.approx(1.0 / 3.0, abs=1e-6)


def test_score_saturates_and_decays_with_excess_mass() -> None:
    # Improving dot quality at fixed mass must never lower the prediction.
    assert calibrated_live_score(40000, 0.10) > calibrated_live_score(40000, 0.05)
    # At fixed w, the optimum is interior: too few dots loses credit, too many
    # pays 0.2 per dot for credit that has already saturated.
    scores = [calibrated_live_score(n, 0.10) for n in range(5_000, 200_001, 5_000)]
    best = int(np.argmax(scores)) * 5_000 + 5_000
    assert 5_000 < best < 200_000
    assert max(scores) > scores[-1]


def test_break_even_rule_is_the_exact_marginal_condition() -> None:
    n, w = 40_000, 0.09
    bar = instrument_break_even(n, w)
    above = calibrated_live_score(n + 1, (n * w + (bar * 1.02)) / (n + 1))
    below = calibrated_live_score(n + 1, (n * w + (bar * 0.98)) / (n + 1))
    base = calibrated_live_score(n, w)
    assert above > base
    assert below < base


def test_predictions_are_bounded_and_monotone_in_w() -> None:
    for n in (1_000, 44_090, 200_000):
        for w in (0.0, 0.02, 0.05, 0.10):
            s = calibrated_live_score(n, w)
            assert 0.0 <= s < 1.0
    assert calibrated_live_score(50_000, 0.0) == 0.0
