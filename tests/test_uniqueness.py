import numpy as np

from gemsdoe40.uniqueness import compare_arrays


def test_identical_raster_is_near_duplicate():
    a = np.zeros((20, 20), dtype=np.float32)
    a[4:8, 5] = 1.0
    footprint = np.ones_like(a, dtype=bool)
    result = compare_arrays(a, a.copy(), footprint, candidate_raw_sha256="a", prior_raw_sha256="a", budget=4)
    assert result["raw_sha256_match"]
    assert result["canonical_pixels_equal_on_footprint"]
    assert result["near_duplicate_by_preregistered_rule"]


def test_disjoint_support_is_not_near_duplicate():
    a = np.zeros((30, 30), dtype=np.float32)
    b = np.zeros_like(a)
    a[2:8, 2] = 1.0
    b[20:26, 27] = 1.0
    result = compare_arrays(a, b, np.ones_like(a, dtype=bool), budget=6)
    assert result["nonzero_support_jaccard"] == 0.0
    assert result["top_budget_containment"] == 0.0
    assert not result["near_duplicate_by_preregistered_rule"]
