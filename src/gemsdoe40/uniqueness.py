"""Hash and correlate a candidate against a pinned inventory of prior TIFFs."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import rasterio

from .raster import canonical_pixel_sha256


NEAR_DUPLICATE_THRESHOLDS = {
    "pearson_abs_gte": 0.995,
    "nonzero_support_jaccard_gte": 0.90,
    "top_budget_containment_gte": 0.90,
}


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _top_positive_indices(values: np.ndarray, k: int) -> np.ndarray:
    flat = np.asarray(values, dtype=np.float32).ravel()
    positive = np.flatnonzero(np.isfinite(flat) & (flat > 0.0))
    if positive.size <= k:
        return positive
    ranks = flat[positive]
    threshold = np.partition(ranks, ranks.size - k)[ranks.size - k]
    higher = positive[ranks > threshold]
    tied = np.sort(positive[ranks == threshold])
    return np.concatenate((higher, tied[: k - higher.size]))


def _correlation(a: np.ndarray, b: np.ndarray) -> float:
    x = np.asarray(a, dtype=np.float64).ravel()
    y = np.asarray(b, dtype=np.float64).ravel()
    if x.size < 2:
        return 0.0
    x0, y0 = x - x.mean(), y - y.mean()
    den = float(np.sqrt(np.dot(x0, x0) * np.dot(y0, y0)))
    if den == 0.0:
        return 1.0 if np.array_equal(x, y) else 0.0
    return float(np.dot(x0, y0) / den)


def compare_arrays(
    candidate: np.ndarray,
    prior: np.ndarray,
    footprint: np.ndarray,
    *,
    candidate_raw_sha256: str | None = None,
    prior_raw_sha256: str | None = None,
    budget: int | None = None,
) -> dict[str, Any]:
    """Compare pixel values and sparse supports on the common template footprint."""
    c = np.asarray(candidate, dtype=np.float32)
    p = np.asarray(prior, dtype=np.float32)
    mask = np.asarray(footprint, dtype=bool)
    if c.shape != p.shape or c.shape != mask.shape:
        raise ValueError("candidate, prior, and footprint shapes must match")
    c = np.nan_to_num(c[mask], nan=0.0, posinf=0.0, neginf=0.0)
    p = np.nan_to_num(p[mask], nan=0.0, posinf=0.0, neginf=0.0)
    corr = _correlation(c, p)
    cm, pm = c > 0.0, p > 0.0
    inter = int(np.count_nonzero(cm & pm))
    union = int(np.count_nonzero(cm | pm))
    jaccard = float(inter / union) if union else 1.0
    containment = float(inter / min(int(cm.sum()), int(pm.sum()))) if cm.any() and pm.any() else 0.0
    k = int(budget or max(int(cm.sum()), int(pm.sum()), 1))
    c_top, p_top = _top_positive_indices(c, k), _top_positive_indices(p, k)
    top_inter = int(np.intersect1d(c_top, p_top, assume_unique=True).size)
    top_containment = float(top_inter / min(c_top.size, p_top.size)) if c_top.size and p_top.size else 0.0
    raw_hash_match = bool(candidate_raw_sha256 and prior_raw_sha256 and candidate_raw_sha256 == prior_raw_sha256)
    pixel_exact = bool(np.array_equal(c, p))
    near = (
        raw_hash_match
        or pixel_exact
        or abs(corr) >= NEAR_DUPLICATE_THRESHOLDS["pearson_abs_gte"]
        or jaccard >= NEAR_DUPLICATE_THRESHOLDS["nonzero_support_jaccard_gte"]
        or top_containment >= NEAR_DUPLICATE_THRESHOLDS["top_budget_containment_gte"]
    )
    return {
        "raw_sha256_match": raw_hash_match,
        "canonical_pixels_equal_on_footprint": pixel_exact,
        "pearson_r_on_footprint": corr,
        "nonzero_intersection": inter,
        "nonzero_union": union,
        "nonzero_support_jaccard": jaccard,
        "nonzero_support_containment": containment,
        "top_budget": k,
        "top_budget_intersection": top_inter,
        "top_budget_containment": top_containment,
        "near_duplicate_by_preregistered_rule": bool(near),
    }


def _artifact_file_name(blob_sha: str) -> str:
    return f"{blob_sha}.tif"


def audit_candidate(
    candidate_path: str | Path,
    template_path: str | Path,
    inventory_path: str | Path,
    prior_cache: str | Path,
) -> dict[str, Any]:
    """Compare against all same-grid unique blobs in the pinned inventory.

    The cache convention is `<git-blob-sha>.tif`. The inventory may contain
    multiple repo paths for the same Git blob; bytes are audited only once.
    """
    candidate_path = Path(candidate_path)
    prior_cache = Path(prior_cache)
    inventory = json.loads(Path(inventory_path).read_text(encoding="utf-8"))
    with rasterio.open(template_path) as template:
        sample = template.read(1)
        footprint = np.isfinite(sample)
        reference = {
            "shape": sample.shape,
            "crs": template.crs,
            "transform": template.transform,
        }
    with rasterio.open(candidate_path) as ds:
        candidate = ds.read(1)
        if candidate.shape != reference["shape"]:
            raise ValueError("candidate does not match template shape")
    candidate_hash = file_sha256(candidate_path)
    candidate_pixel_hash = canonical_pixel_sha256(candidate, ~footprint)
    results: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    for entry in inventory.get("unique_rasters", []):
        blob = entry["git_blob_sha"]
        cached = prior_cache / _artifact_file_name(blob)
        if not cached.exists():
            skipped.append({"git_blob_sha": blob, "reason": "not present in prior cache"})
            continue
        with rasterio.open(cached) as prior_ds:
            same_grid = (
                prior_ds.count == 1
                and prior_ds.height == reference["shape"][0]
                and prior_ds.width == reference["shape"][1]
                and prior_ds.crs == reference["crs"]
                and prior_ds.transform == reference["transform"]
            )
            if not same_grid:
                skipped.append({"git_blob_sha": blob, "reason": "not a single-band exact-template grid"})
                continue
            prior = prior_ds.read(1)
        prior_hash = file_sha256(cached)
        expected_hash = entry.get("sha256")
        if expected_hash and prior_hash != expected_hash:
            raise ValueError(f"prior raster SHA-256 mismatch for {blob}: {prior_hash} != {expected_hash}")
        metrics = compare_arrays(
            candidate,
            prior,
            footprint,
            candidate_raw_sha256=candidate_hash,
            prior_raw_sha256=prior_hash,
            budget=int(np.count_nonzero(candidate[footprint])),
        )
        results.append({
            "git_blob_sha": blob,
            "sha256": prior_hash,
            "bytes": cached.stat().st_size,
            "paths": entry.get("artifacts", []),
            **metrics,
        })
    results.sort(key=lambda row: (row["near_duplicate_by_preregistered_rule"], abs(row["pearson_r_on_footprint"]), row["nonzero_support_jaccard"]), reverse=True)
    closest = results[:10]
    return {
        "candidate": {
            "path": str(candidate_path),
            "sha256": candidate_hash,
            "canonical_pixels_sha256": candidate_pixel_hash,
            "nonzero_in_footprint": int(np.count_nonzero(candidate[footprint])),
        },
        "thresholds": NEAR_DUPLICATE_THRESHOLDS,
        "inventory_unique_blobs": len(inventory.get("unique_rasters", [])),
        "same_grid_comparisons": len(results),
        "different_grid_or_band_count": len(skipped),
        "near_duplicate_count": sum(bool(row["near_duplicate_by_preregistered_rule"]) for row in results),
        "uniqueness_pass": not any(row["near_duplicate_by_preregistered_rule"] for row in results),
        "closest_prior_outputs": closest,
        "comparisons": results,
        "skipped": skipped,
    }
