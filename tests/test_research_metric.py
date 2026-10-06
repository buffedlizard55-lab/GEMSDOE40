import numpy as np

from gemsdoe40.research_metric import metric_components, spatial_block_components


def test_exact_hit_scores_one():
    truth = np.zeros((9, 9), dtype=bool)
    truth[4, 4] = True
    pred = truth.astype(np.float32)
    score = metric_components(pred, truth)
    assert score.tp == 1.0
    assert score.fp == 0.0
    assert score.fn == 0.0
    assert np.isclose(score.score, 1.0)


def test_one_pixel_near_miss_uses_triangular_weights():
    truth = np.zeros((9, 9), dtype=bool)
    truth[4, 4] = True
    pred = np.zeros((9, 9), dtype=np.float32)
    pred[4, 5] = 1.0
    score = metric_components(pred, truth, radius_px=3.0)
    assert np.isclose(score.tp, 2.0 / 3.0)
    assert np.isclose(score.fp, 1.0 / 3.0)
    assert np.isclose(score.fn, 1.0 / 3.0)
    assert np.isclose(score.score, 2.0 / 3.0)


def test_predictions_past_three_pixels_are_full_false_positive():
    truth = np.zeros((12, 12), dtype=bool)
    truth[2, 2] = True
    pred = np.zeros_like(truth, dtype=np.float32)
    pred[2, 7] = 1.0
    score = metric_components(pred, truth)
    assert score.tp == 0.0
    assert score.fp == 1.0
    assert score.fn == 1.0


def test_exact_known_pixel_mask_and_spatial_blocks():
    truth = np.zeros((24, 24), dtype=bool)
    truth[5, 5] = True
    truth[18, 18] = True
    valid = np.ones_like(truth)
    pred = truth.astype(np.float32)
    score, blocks = spatial_block_components(pred, truth, valid, n_rows=2, n_cols=2, guard=1)
    assert np.isclose(score.score, 1.0)
    assert len(blocks) == 4
    assert sum(block["n_truth"] for block in blocks) == 2
