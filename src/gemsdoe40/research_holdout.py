"""Frozen SGMC-proxy holdout loader, scorer, and incumbent comparison."""
from __future__ import annotations

from pathlib import Path
from typing import Any
import hashlib

import numpy as np
import rasterio

from .research_metric import spatial_block_components
from gems40.pins import PINNED_FILES

# Historical H2-B evidence and the current proxy are distinct file versions. Keep the default
# historical pin for old runner compatibility; current experiments must pass the current pin
# explicitly so the two instruments can never be silently conflated.
EXPECTED_PROXY_SHA256 = "26d142c4c93282cd94f6950ab96f22aeff59fbbea523d43d662e76fa1b161b5c"
CURRENT_PROXY_SHA256 = PINNED_FILES["derived_sgmc_faults_100m_u8.tif"]["sha256"]


def read_proxy_truth(
    proxy_path: str | Path,
    template_path: str | Path,
    labels_path: str | Path,
    *,
    expected_proxy_sha256: str = EXPECTED_PROXY_SHA256,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """Build an exact-pixel off-catalogue mask; do not buffer known faults.

    The default preserves the historical H2-B pin. Current experiments must explicitly pass
    ``CURRENT_PROXY_SHA256`` rather than inheriting the historical instrument by accident.
    """
    proxy_sha256 = hashlib.sha256(Path(proxy_path).read_bytes()).hexdigest()
    if proxy_sha256 != expected_proxy_sha256:
        raise ValueError(f"SGMC proxy hash mismatch: {proxy_sha256} != pinned {expected_proxy_sha256}")
    with rasterio.open(template_path) as template, rasterio.open(labels_path) as labels_ds, rasterio.open(proxy_path) as proxy_ds:
        sample = template.read(1)
        valid = np.isfinite(sample)
        labels = labels_ds.read(1)
        proxy = proxy_ds.read(1)
        for ds, name in ((labels_ds, "labels"), (proxy_ds, "proxy")):
            if ds.shape != template.shape or ds.crs != template.crs or ds.transform != template.transform:
                raise ValueError(f"{name} raster is not exactly aligned to the sample template")
        if proxy_ds.count != 1 or labels_ds.count != 1:
            raise ValueError("proxy and known labels must be single-band rasters")
        if not np.isin(np.unique(proxy[valid]), [0, 1]).all():
            raise ValueError("proxy has values other than 0/1 inside the template footprint")
        known = (labels == 1) & valid
        truth = (proxy == 1) & valid & (labels != 1)
        info = {
            "proxy_status": "owner-derived SGMC mirror; proxy instrument only, not organizer truth",
            "proxy_sha256": proxy_sha256,
            "proxy_generation": "historical_h2b_26d142" if proxy_sha256 == EXPECTED_PROXY_SHA256 else ("current_643cbe" if proxy_sha256 == CURRENT_PROXY_SHA256 else "explicitly_pinned_other"),
            "proxy_dtype": proxy_ds.dtypes[0],
            "proxy_nodata": proxy_ds.nodata,
            "template_valid_pixels": int(valid.sum()),
            "template_outside_pixels": int((~valid).sum()),
            "known_label_pixels_exact": int(known.sum()),
            "sgmc_positive_pixels_in_footprint": int(np.count_nonzero((proxy == 1) & valid)),
            "sgmc_positive_overlap_exact_known_pixels": int(np.count_nonzero((proxy == 1) & known)),
            "off_catalogue_truth_pixels_after_exact_mask": int(truth.sum()),
            "mask_rule": "exclude only labels == 1 at those exact pixels; no distance buffer",
            "holdout_blocks": "4 rows x 6 columns, contiguous; 3-cell guard at internal boundaries",
            "official_metric": "published 300 m triangular distance-weighted Tversky; alpha=0.2, beta=0.8",
        }
        return truth, valid, labels, info


def score_array_on_proxy(
    prediction: np.ndarray,
    truth: np.ndarray,
    valid: np.ndarray,
    labels: np.ndarray,
) -> dict[str, Any]:
    """Score one prediction on already-loaded frozen proxy arrays."""
    pred = np.asarray(prediction, dtype=np.float32).copy()
    if pred.shape != truth.shape:
        raise ValueError("prediction shape does not match proxy grid")
    # DrivenData's known mask is pixel-exact. Apply no buffer around it.
    pred[(labels == 1) & valid] = 0.0
    pooled, blocks = spatial_block_components(
        pred,
        truth,
        valid,
        n_rows=4,
        n_cols=6,
        guard=3,
        alpha=0.2,
        beta=0.8,
        radius_px=3.0,
    )
    return {"pooled": pooled.to_dict(), "blocks": blocks}


def score_candidate_on_proxy(
    prediction: np.ndarray,
    proxy_path: str | Path,
    template_path: str | Path,
    labels_path: str | Path,
    *,
    expected_proxy_sha256: str = EXPECTED_PROXY_SHA256,
) -> tuple[dict[str, Any], np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """Load and score a prediction on the locked 24-block proxy holdout."""
    truth, valid, labels, info = read_proxy_truth(
        proxy_path,
        template_path,
        labels_path,
        expected_proxy_sha256=expected_proxy_sha256,
    )
    score = score_array_on_proxy(prediction, truth, valid, labels)
    return score, truth, valid, labels, info


def promotion_gate(candidate: dict[str, Any], incumbent: dict[str, Any]) -> dict[str, Any]:
    """Apply the registered no-slot-unless-better rule to 24 guarded blocks."""
    cblocks = candidate["blocks"]
    iblocks = incumbent["blocks"]
    if len(cblocks) != 24 or len(iblocks) != 24:
        raise ValueError("promotion gate requires exactly 24 spatial blocks")
    deltas = [float(c["score"] - i["score"]) for c, i in zip(cblocks, iblocks)]
    wins = sum(delta > 0.0 for delta in deltas)
    pooled_delta = float(candidate["pooled"]["score"] - incumbent["pooled"]["score"])
    passed = pooled_delta > 0.0 and wins >= 18 and min(deltas) >= -0.005
    return {
        "passed": bool(passed),
        "decision": "ELIGIBLE_FOR_REVIEW; still requires human portal confirmation" if passed else "HOLD; DO NOT SUBMIT",
        "candidate_pooled_dti": candidate["pooled"]["score"],
        "incumbent_pooled_dti": incumbent["pooled"]["score"],
        "pooled_delta": pooled_delta,
        "blocks_won": int(wins),
        "blocks_required": 18,
        "minimum_block_delta": float(min(deltas)),
        "maximum_block_delta": float(max(deltas)),
        "registered_minimum_block_delta": -0.005,
        "per_block_deltas": deltas,
    }
