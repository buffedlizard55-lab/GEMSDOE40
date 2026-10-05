import numpy as np
import pytest

from gemsdoe40.euler import emit_all_supported, emit_top_budget, solve_euler_patch


def test_si0_euler_patch_recovers_constructed_source():
    n = 10
    res = 100.0
    offsets = (np.arange(n) + 0.5 - n / 2.0) * res
    x = np.broadcast_to(offsets[None, :], (n, n))
    y = np.broadcast_to((-offsets)[:, None], (n, n))
    x0, y0, z0 = 120.0, -60.0, 700.0
    rng = np.random.default_rng(51)
    tx = rng.normal(size=(n, n))
    ty = rng.normal(size=(n, n))
    tz = ((x - x0) * tx + (y - y0) * ty) / z0
    got_x, got_y, got_z, residual, condition = solve_euler_patch(tx, ty, tz)
    assert np.allclose([got_x, got_y, got_z], [x0, y0, z0], rtol=1e-9, atol=1e-7)
    assert residual < 1e-7
    assert np.isfinite(condition)


def test_emitter_is_binary_sparse_and_uses_row_major_ties():
    kde = np.ones((4, 5), dtype=np.float32)
    valid = np.ones_like(kde, dtype=bool)
    known = np.zeros_like(kde, dtype=np.int8)
    pred, info = emit_top_budget(kde, valid, known, budget=5)
    assert np.count_nonzero(pred) == 5
    assert np.array_equal(np.flatnonzero(pred), np.arange(5))
    assert set(np.unique(pred)) == {0.0, 1.0}
    assert info["emitted_pixels"] == 5


def test_all_support_emitter_keeps_only_positive_non_catalogue_kde():
    kde = np.array([[1.0, 0.2, 0.0], [0.4, 0.3, 0.1]], dtype=np.float32)
    valid = np.ones_like(kde, dtype=bool)
    known = np.zeros_like(kde, dtype=np.int8)
    known[0, 1] = 1
    pred, info = emit_all_supported(kde, valid, known)
    assert np.count_nonzero(pred) == 4
    assert pred[0, 1] == 0.0
    assert pred[0, 0] == 1.0
    assert info["emitted_pixels"] == 4


def test_fixed_budget_emitter_refuses_to_backfill_when_support_is_short():
    kde = np.array([[1.0, 0.0], [0.0, 0.5]], dtype=np.float32)
    valid = np.ones_like(kde, dtype=bool)
    known = np.zeros_like(kde, dtype=np.int8)
    with pytest.raises(RuntimeError, match="refuse to backfill"):
        emit_top_budget(kde, valid, known, budget=3)
