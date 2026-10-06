"""Analytic tests for the H8 lineament-aware depth-clustering construction.

Everything here is synthetic: no competition raster, no label, no prior
prediction.  The tests pin (a) the geometry convention (a synthetic straight
lineament must be recovered as a lineament *along* its own strike, not across
it), (b) the depth-consensus rule (cross-window agreement is required before a
solution keeps weight), (c) the maximum-separation emission (no two dots closer
than the spacing, and the highest-density cell is always kept first), and
(d) the exact distance-weighted Tversky identity used by the audit.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

SPEC = importlib.util.spec_from_file_location("run_h8_lineament", ROOT / "scripts/run_h8_lineament.py")
h8 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(h8)

from gemsdoe40.contact_euler import CLOUD_DTYPE  # noqa: E402


def synthetic_cloud(rows, cols, depths, windows):
    cloud = np.zeros(len(rows), dtype=h8.CLOUD_EXT_DTYPE)
    cloud["row"] = rows
    cloud["col"] = cols
    cloud["depth_m"] = depths
    cloud["depth_se_m"] = 25.0
    cloud["residual"] = 0.05
    cloud["condition"] = 2.0
    cloud["window"] = np.asarray(windows, dtype=np.uint8)
    return cloud


def test_local_axes_recovers_strike_of_a_straight_train():
    """A NE-trending train (row and col both increasing) must come back as 45°."""
    steps = np.arange(40)
    rows = 500.0 + steps * 1.0
    cols = 800.0 + steps * 1.0
    cloud = synthetic_cloud(rows, cols, np.full(steps.size, 600.0), np.full(steps.size, 15))
    direction, coherence = h8.local_axes(cloud, cell_m=100.0)
    middle = 20
    assert abs(np.rad2deg(direction[middle]) - 45.0) < 5.0, direction[middle]
    assert coherence[middle] > 0.9


def test_anisotropic_kernel_smears_along_the_train_not_across_it():
    """The deposited density must be anisotropic in the direction of the train."""
    steps = np.arange(41)
    rows = 600.0 + steps * 1.0          # 45 degrees in (row, col)
    cols = 600.0 + steps * 1.0
    cloud = synthetic_cloud(rows, cols, np.full(steps.size, 500.0), np.full(steps.size, 15))
    cloud["lineament_coherence"] = 1.0
    row0 = int(rows[0])
    col0 = int(cols[0])
    field = h8.anisotropic_kde(cloud["row"], cloud["col"], np.ones(steps.size),
                               np.full(steps.size, np.pi / 4), cloud["lineament_coherence"],
                               (2000, 2000))
    # compare mass deposited along the train direction with the mass across it
    along = sum(field[row0 + k, col0 + k] for k in range(-6, 7))
    across = sum(field[row0 + k, col0 - k] for k in range(-6, 7))
    assert along > 5.0 * across, (along, across)


def test_consensus_weights_require_cross_window_agreement():
    """Identical xy but a single window -> no weight; two windows -> weight."""
    rows = np.full(6, 400.0) + np.arange(6)
    cols = np.full(6, 700.0) + np.arange(6)
    depths = np.full(6, 800.0)
    single = synthetic_cloud(rows, cols, depths, np.full(6, 15))
    assert np.all(h8.consensus_weights(single, cell_m=100.0) == 0.0)
    windows = np.array([15, 15, 15, 21, 21, 21], dtype=np.uint8)
    mixed = synthetic_cloud(rows, cols, depths, windows)
    weight = h8.consensus_weights(mixed, cell_m=100.0)
    assert weight.max() > 0.0 and weight.min() >= 0.0
    # shallower solutions must outrank deeper ones at equal support
    deep = synthetic_cloud(rows, cols, depths * 4.0, windows)
    assert h8.consensus_weights(deep, cell_m=100.0).max() < weight.max()


def test_value_ranked_nms_respects_spacing_and_priority():
    field = np.zeros((40, 40), dtype=np.float64)
    field[5, 5] = 1.0
    field[5, 7] = 0.9        # 2 px away -> suppressed
    field[20, 20] = 0.8      # far away -> kept
    eligible = np.zeros_like(field, dtype=bool)
    eligible[field > 0] = True
    mask = h8.value_ranked_nms(field, eligible, spacing_px=2.8, budget=10)
    assert mask[5, 5] and mask[20, 20] and not mask[5, 7]
    assert mask.sum() == 2


def test_audit_metric_identity_matches_the_official_algebra():
    """The audit's surrogate must equal the metric form in gemsdoe40.emission."""
    from gemsdoe40.emission import credit_of_binary, pred_credit_of_binary
    spec = importlib.util.spec_from_file_location("audit_h8_lineament",
                                                 ROOT / "scripts/audit_h8_lineament.py")
    audit = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit)

    rng = np.random.default_rng(3)
    truth = np.zeros((40, 40), dtype=bool)
    truth[10:14, 10:25] = True
    prediction = np.zeros((40, 40), dtype=np.float32)
    prediction[9:13, 12:20] = 1.0
    prediction[30, 30] = 1.0

    from scipy import ndimage
    distance = ndimage.distance_transform_edt(~truth)
    kernel = audit.kernel_credit(distance)
    result = audit.stats(prediction, truth, kernel=kernel)
    tp = credit_of_binary(prediction, truth)
    pred_credit = pred_credit_of_binary(prediction, truth)
    assert result["tp"] == pytest.approx(tp, rel=1e-6)
    assert result["pred_credit"] == pytest.approx(pred_credit, rel=1e-6)
    fn = truth.sum() - tp
    fp = prediction.sum() - pred_credit
    assert result["dti_surrogate"] == pytest.approx(tp / (tp + 0.2 * fp + 0.8 * fn), rel=1e-6)
