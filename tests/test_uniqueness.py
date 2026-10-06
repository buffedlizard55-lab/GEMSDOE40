from __future__ import annotations

import numpy as np

from gemsdoe40.uniqueness import pearson, jaccard_positive, NEAR_DUPLICATE_CORR


def test_pearson_identical_is_one():
    a = np.arange(12, dtype=float).reshape(3, 4)
    assert pearson(a, a) == np.float64(1.0) or abs(pearson(a, a) - 1) < 1e-12


def test_pearson_orthogonal_low():
    a = np.zeros((8, 8))
    b = np.zeros((8, 8))
    a[0, 0] = 1
    b[7, 7] = 1
    assert abs(pearson(a, b)) < 0.3


def test_jaccard_disjoint_zero():
    a = np.zeros((5, 5))
    b = np.zeros((5, 5))
    a[0, 0] = 1
    b[4, 4] = 1
    assert jaccard_positive(a, b) == 0.0


def test_thresholds_are_strict():
    assert NEAR_DUPLICATE_CORR == 0.85
