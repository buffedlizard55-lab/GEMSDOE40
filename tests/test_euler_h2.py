import numpy as np
import pytest

from gemsdoe40.euler import EulerCloud, persistence_cloud, upward_continue_with_vertical


def _cloud(row, col, depth, residual=0.1):
    row = np.asarray(row, dtype=np.float32)
    return EulerCloud(
        row,
        np.asarray(col, dtype=np.float32),
        np.asarray(depth, dtype=np.float32),
        np.full(row.size, residual, dtype=np.float32),
        np.ones(row.size, dtype=np.float32),
    )


def test_persistence_requires_three_heights_and_returns_weighted_base_points():
    clouds = {
        0: _cloud([100.0], [200.0], [1500.0]),
        500: _cloud([100.5], [200.0], [1510.0]),
        1000: _cloud([101.0], [200.5], [1490.0]),
        2000: _cloud([], [], []),
    }
    pairs, hits = persistence_cloud(clouds)
    assert len(pairs.row) == 1
    assert hits.tolist() == [3]
    assert pairs.row[0] == 100.0
    assert pairs.col[0] == 200.0
    terms = [
        (0.1 / 0.20) ** 2,
        (0.1 / 0.20) ** 2 + (50.0 / 300.0) ** 2 + (10.0 / 1000.0) ** 2,
        (0.1 / 0.20) ** 2 + ((12500.0**0.5) / 300.0) ** 2 + (10.0 / 1000.0) ** 2,
    ]
    expected_weight = (3 / 4) * np.exp(-0.5 * np.mean(terms))
    assert pairs.weight[0] == pytest.approx(expected_weight, rel=1e-6)


def test_upward_continuation_preserves_constant_and_has_zero_vertical_derivative():
    field = np.full((24, 32), 7.5, dtype=np.float32)
    continued, vertical = upward_continue_with_vertical(field, 1000.0, pad_cells=8)
    assert np.allclose(continued, field, atol=1e-5)
    assert np.max(np.abs(vertical)) < 1e-5


def test_persistence_rejects_two_height_only_solution():
    clouds = {
        0: _cloud([100.0], [200.0], [1500.0]),
        500: _cloud([100.5], [200.0], [1510.0]),
        1000: _cloud([], [], []),
        2000: _cloud([], [], []),
    }
    pairs, hits = persistence_cloud(clouds)
    assert len(pairs.row) == 0
    assert hits.size == 0
