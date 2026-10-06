"""Cross-checks between the fast measurement layer and the pinned metric."""

from __future__ import annotations

import numpy as np
import pytest

from gemsdoe40 import measure
from gemsdoe40.research_metric import metric_components, spatial_block_components


def _case(seed: int = 0, shape=(70, 90), mass: int = 500):
    rng = np.random.default_rng(seed)
    truth = rng.random(shape) < 0.02
    pred = np.zeros(shape, dtype=np.float32)
    pred.ravel()[rng.choice(np.arange(truth.size), size=mass, replace=False)] = 1.0
    return pred, truth


@pytest.mark.parametrize("seed", [0, 1, 2, 3])
def test_components_match_research_metric(seed: int) -> None:
    pred, truth = _case(seed)
    fast = measure.credit_components(pred, truth)
    slow = metric_components(pred, truth)
    assert fast["tp"] == pytest.approx(slow.tp, rel=1e-6)
    assert fast["fp"] == pytest.approx(slow.fp, rel=1e-6)
    assert fast["fn"] == pytest.approx(slow.fn, rel=1e-6)
    assert fast["dti"] == pytest.approx(slow.score, rel=1e-6)


@pytest.mark.parametrize("seed", [4, 5])
def test_components_match_with_mask(seed: int) -> None:
    pred, truth = _case(seed)
    rng = np.random.default_rng(seed + 100)
    mask = rng.random(truth.shape) < 0.01
    fast = measure.credit_components(pred, truth, mask=mask)
    slow = metric_components(pred, truth, valid=~mask)
    assert fast["dti"] == pytest.approx(slow.score, rel=1e-6)


def test_blocked_matches_research_metric_blocks() -> None:
    pred, truth = _case(11, shape=(80, 120), mass=800)
    valid = np.ones(pred.shape, dtype=bool)
    fast = measure.blocked_components(pred, truth, valid=valid, n_rows=4, n_cols=6, guard=3)
    slow_pooled, slow_blocks = spatial_block_components(pred, truth, valid, n_rows=4, n_cols=6, guard=3)
    assert fast["pooled"]["dti"] == pytest.approx(slow_pooled.score, rel=1e-6)
    for a, b in zip(fast["blocks"], slow_blocks):
        assert (a["row"], a["col"]) == (b["row"], b["col"])
        assert a["dti"] == pytest.approx(b["score"], rel=1e-5, abs=1e-9)


def test_kernel_field_is_the_credit_each_truth_pixel_receives() -> None:
    pred, truth = _case(21)
    K = measure.kernel_field(pred)
    rng = np.random.default_rng(0)
    subset = truth & (rng.random(truth.shape) < 0.5)
    assert K[subset].sum() == pytest.approx(measure.credit_components(pred, subset)["tp"], rel=1e-6)


def test_same_mass_random_matches_mass() -> None:
    pred, _ = _case(2)
    valid = np.ones(pred.shape, dtype=bool)
    rnd = measure.same_mass_random(pred, valid, np.random.default_rng(1))
    assert int(rnd.sum()) == int((pred > 0).sum())
