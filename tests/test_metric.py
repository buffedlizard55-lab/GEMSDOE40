"""Reproduce the official worked example and the algebraic identities."""
from __future__ import annotations

import numpy as np
import pytest

from gemsdoe40.metric import dti_from_components, dti_exact, dti_binary, kernel, marginal_bar


def test_official_worked_example_components():
    # Page 967: TP_w=3.00, FP_w=1.89, FN_w=2.00 → TI_w = 0.60
    dti = dti_from_components(3.00, 1.89, 2.00)
    assert dti == pytest.approx(3.00 / (3.00 + 0.2 * 1.89 + 0.8 * 2.00), rel=1e-12)
    assert dti == pytest.approx(0.60, abs=5e-3)


def test_kernel_at_lattice_offsets():
    assert kernel(0) == 1.0
    assert kernel(3) == 0.0
    assert kernel(1.5) == pytest.approx(0.5)
    assert kernel(4) == 0.0


def test_perfect_overlap_binary():
    g = np.zeros((9, 9), dtype=bool)
    g[4, 4] = True
    p = g.astype(float)
    res = dti_binary(g, g)
    assert res["dti"] == pytest.approx(1.0)
    assert res["tp"] == pytest.approx(1.0)
    assert res["fp"] == pytest.approx(0.0)
    assert res["fn"] == pytest.approx(0.0)


def test_soft_and_binary_agree_on_01():
    rng = np.random.default_rng(0)
    p = (rng.random((12, 12)) > 0.8).astype(float)
    g = rng.random((12, 12)) > 0.8
    a = dti_exact(p, g)
    b = dti_binary(p > 0, g)
    assert a["dti"] == pytest.approx(b["dti"], rel=1e-10)
    assert a["tp"] == pytest.approx(b["tp"], rel=1e-10)


def test_marginal_bar_is_alpha_times_dti():
    assert marginal_bar(0.2708) == pytest.approx(0.05416)
    assert marginal_bar(0.2778) == pytest.approx(0.05556)
