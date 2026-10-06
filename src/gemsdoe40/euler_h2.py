"""Preregistered H2 wrapper: multi-height persistence of TMI Euler solutions."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import rasterio

from .research_euler import (
    EMISSION_BUDGET,
    kde_from_pairs,
    persistence_cloud,
    solve_euler_cloud,
    emit_all_supported,
    upward_continue_with_vertical,
    _gradient_components,
)
from .raster import band_index_by_name

HEIGHTS_M = (0, 500, 1000, 2000)


def run_h2(
    features_path: str | Path,
    template_path: str | Path,
    labels_path: str | Path,
    *,
    support_only: bool = False,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Generate H2/H2-B without using holdout labels or leaderboard scores."""
    features_path, template_path, labels_path = Path(features_path), Path(template_path), Path(labels_path)
    with rasterio.open(template_path) as sample, rasterio.open(labels_path) as labels_ds, rasterio.open(features_path) as features:
        template = sample.read(1)
        footprint = np.isfinite(template)
        labels = labels_ds.read(1)
        if labels.shape != template.shape or features.shape != template.shape:
            raise ValueError("sample, labels, and features must have identical shape")
        band = band_index_by_name(features, "tmi")
        tmi = features.read(band).astype(np.float32, copy=False)
        transform = features.transform
        if sample.crs != features.crs or sample.transform != features.transform:
            raise ValueError("sample and feature georeferencing differ")
    input_valid = footprint & np.isfinite(tmi) & (np.abs(tmi) < 1e30)

    clouds = {}
    for height in HEIGHTS_M:
        continued, tz = upward_continue_with_vertical(tmi, float(height))
        tx, ty = _gradient_components(continued)
        clouds[height] = solve_euler_cloud(tx, ty, tz, input_valid, height_m=float(height))
    persistent, hits = persistence_cloud(clouds)
    kde = kde_from_pairs(persistent, footprint.shape)
    if support_only:
        pred, emission = emit_all_supported(kde, footprint, labels)
        hypothesis_id = "H2-B"
    else:
        from .research_euler import emit_top_budget  # keep H1/H2 on the same preregistered emitter
        pred, emission = emit_top_budget(kde, footprint, labels, budget=EMISSION_BUDGET)
        hypothesis_id = "H2"
    pred[labels == 1] = 0.0
    pred[~footprint] = 0.0
    summary: dict[str, Any] = {
        "hypothesis_id": hypothesis_id,
        "transform": list(transform)[:6],
        "shape": list(footprint.shape),
        "feature_band": {"band": band, "band_name": "tmi"},
        "upward_continuation_heights_m": list(HEIGHTS_M),
        "vertical_derivative": "padded reflect +|k| derivative at every upward-continuation height; positive down",
        "common_real_input_pixels": int(input_valid.sum()),
        "exact_known_label_pixels_suppressed": int(np.count_nonzero((labels == 1) & footprint)),
        "clouds_by_height": {str(h): cloud.to_summary() for h, cloud in clouds.items()},
        "persistent_solutions": {
            **persistent.to_summary(),
            "matches_at_3_heights": int(np.count_nonzero(hits == 3)),
            "matches_at_4_heights": int(np.count_nonzero(hits == 4)),
        },
        "emission": emission,
    }
    return pred, summary
