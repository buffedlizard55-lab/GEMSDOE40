"""Fail-closed H4 raw-output novelty checks; no random pixel sampling."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import rasterio
from rasterio.warp import Resampling, reproject

from .raster import canonical_pixel_sha256, validate_candidate
from .research_uniqueness import file_sha256

RULE = {"pearson_abs_gte": .85, "top_budget": 37654,
        "top_jaccard_gte": .50, "top_containment_gte": .80}

# Forensically inspected exact historical bytes, NOT a blanket exemption for
# unreadable files, small rasters, or filenames. See the public IFD/payload
# record docs/data/prior-nonprediction-forensics-20261006.json.
SPATIAL_FORMAT_FIXTURE_SHA256 = "60543f8ea9d4ee432b018832337c3cf1fe0db4f71b87d15fea71869fe4f23638"
NONPREDICTION_SHA256 = {
    "de5e3531caedbaed73735d7de4f426f86f79c42c7f86aae5edacdbe24e5158b6":
        "110-byte TIFF header ends immediately after its IFD; no StripOffsets, byte counts or pixel payload. Not a raw prediction raster.",
    "bdf5dbb6f1f8ef758dfdf4c4fd3e086845d30720af91e2e4c6600d759396039e":
        "138-byte ungeoreferenced 1x1 demo; exactly one float32 value 0.7310000061988831, not a GeoDAWN prediction.",
}


# Two historical "nan" files incorrectly mask many INSIDE-footprint zero
# cells. Their separately published all-finite companions contain exactly the
# same entire array with those NaNs represented by zeros. Every companion is
# itself hashed and fully compared in the corpus. No unknown value is guessed.
NAN_STORAGE_COUNTERPARTS = {
    "fbef100fbadc57ec1c663d523101114ba312d81518ffd8fd21229d91dacf3835":
        ("3ec3182c9ebc01eeca86672de1559d2bdda66162", "1f729e73c8ff1d23b1bafe0bc6474198294863a40efd025a25152d4b22cf2325"),
    "29f5c1b4edf4832408f4a074cbca91b5965c5aba4d83c911362e59a6cd9e8e68":
        ("ca1434487c5f9a0fcae74b75ed8cedfad3631bf3", "10e2e7a195dd8109eff80bbcc55a944476982460427ef740c0c90d82d507c573"),
}


def nan_zero_storage_equivalent(partial: np.ndarray, complete: np.ndarray) -> bool:
    """Prove a published storage relationship; never mutate comparison values."""
    return bool(partial.shape == complete.shape and np.isfinite(complete).all()
                and not np.isinf(partial).any()
                and np.array_equal(np.where(np.isnan(partial), 0, partial), complete))


def top_indices(values: np.ndarray, k: int) -> np.ndarray:
    """Positive top-k, with deterministic ascending-index tie break."""
    v = np.asarray(values).ravel()
    ids = np.flatnonzero(np.isfinite(v) & (v > 0))
    if len(ids) <= k:
        return ids
    q = v[ids]
    threshold = np.partition(q, len(q) - k)[len(q) - k]
    higher = ids[q > threshold]
    ties = ids[q == threshold]
    return np.concatenate((higher, ties[:k-len(higher)]))


def correlation(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.size < 2:
        return 0.0
    equal = np.array_equal(x, y)
    x = x - x.mean(); y = y - y.mean()
    denom = np.linalg.norm(x) * np.linalg.norm(y)
    return float(np.dot(x, y) / denom) if denom > 0 else (1.0 if equal else 0.0)


def compare_vectors(candidate: np.ndarray, prior: np.ndarray, *, budget: int = RULE["top_budget"]) -> dict[str, Any]:
    c, p = (np.asarray(a, dtype=np.float32).ravel() for a in (candidate, prior))
    if c.shape != p.shape or not c.size or not (np.isfinite(c).all() and np.isfinite(p).all()):
        raise ValueError("nonempty equal-length finite vectors required; do not silently sanitize raw outputs")
    r = correlation(c, p)
    ct, pt = top_indices(c, budget), top_indices(p, budget)
    inter = len(np.intersect1d(ct, pt, assume_unique=True))
    union = len(ct) + len(pt) - inter
    jac = inter / union if union else 0.0
    contain = inter / min(len(ct), len(pt)) if len(ct) and len(pt) else 0.0
    cm, pm = c > 0, p > 0
    support_inter, support_union = int(np.count_nonzero(cm & pm)), int(np.count_nonzero(cm | pm))
    support_j = support_inter / support_union if support_union else 0.0
    # Reproduce the historic rule separately; broad positive support is not our
    # definition of continuous-raster duplication.
    old_k = max(int(cm.sum()), 1)
    old_ct, old_pt = top_indices(c, old_k), top_indices(p, old_k)
    old_inter = len(np.intersect1d(old_ct, old_pt, assume_unique=True))
    old_contain = old_inter / min(len(old_ct), len(old_pt)) if len(old_ct) and len(old_pt) else 0.0
    exact = bool(np.array_equal(c, p))
    hot = {}
    for q in (.95, .99):
        ca = c >= np.quantile(c[cm], q) if cm.any() else np.zeros(c.shape, bool)
        pa = p >= np.quantile(p[pm], q) if pm.any() else np.zeros(p.shape, bool)
        u = np.count_nonzero(ca | pa)
        hot[str(q)] = float(np.count_nonzero(ca & pa) / u) if u else 0.0
    return {
        "pixels_compared": int(c.size), "pearson_r": r, "exact_values_equal": exact,
        "top_budget": budget, "candidate_top_pixels": len(ct), "prior_top_pixels": len(pt),
        "top_intersection": inter, "top_jaccard": float(jac), "top_containment": float(contain),
        "nonzero_support_jaccard": float(support_j), "hot_support_jaccards": hot,
        "near_duplicate_h4": bool(exact or abs(r) >= RULE["pearson_abs_gte"] or jac >= RULE["top_jaccard_gte"] or contain >= RULE["top_containment_gte"]),
        "historic_h2b_rule_near_duplicate": bool(exact or abs(r) >= .995 or support_j >= .9 or old_contain >= .9),
        "historic_top_budget_containment": float(old_contain),
    }


def audit(candidate_path: Path, sample_path: Path, inventory_path: Path, cache: Path) -> dict:
    format_receipt = validate_candidate(candidate_path, sample_path)
    if not format_receipt["valid"]:
        raise ValueError(f"candidate must pass exact format before novelty: {format_receipt['issues']}")
    inventory = json.loads(inventory_path.read_text())
    entries = inventory.get("unique_rasters", [])
    if not entries or inventory.get("errors"):
        raise ValueError("empty/incomplete inventory cannot establish uniqueness")
    with rasterio.open(sample_path) as s:
        foot = np.isfinite(s.read(1)); shape, crs, transform = s.shape, s.crs, s.transform
    with rasterio.open(candidate_path) as s:
        candidate = s.read(1)
    results, issues = [], []
    sha = file_sha256(candidate_path)
    canonical = canonical_pixel_sha256(candidate, ~foot)
    for i, e in enumerate(entries, 1):
        path = cache / f"{e['git_blob_sha']}.tif"
        if not path.is_file():
            issues.append({"blob": e["git_blob_sha"], "issue": "missing prior bytes"}); continue
        raw_sha = file_sha256(path)
        data = path.read_bytes()
        actual_blob = hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()
        if len(data) != e["bytes"] or raw_sha != e["sha256"] or actual_blob != e["git_blob_sha"]:
            raise ValueError(f"corrupt prior bytes: {path}")
        del data
        row = {"git_blob_sha": e["git_blob_sha"], "sha256": raw_sha, "bytes": e["bytes"],
               "artifact": e["artifacts"][0], "file_hash_equal": raw_sha == sha}
        if raw_sha in NONPREDICTION_SHA256:
            row.update(exact_grid=False, near_duplicate_h4=False,
                       nonprediction_artifact=True, correlation_applicable=False,
                       comparison_mode="exact-byte forensic classification, NOT a raw-map correlation",
                       classification=NONPREDICTION_SHA256[raw_sha])
            results.append(row)
            continue
        try:
            dataset = rasterio.open(path)
        except rasterio.errors.RasterioIOError as exc:
            row.update(exact_grid=False, near_duplicate_h4=False, unreadable=True)
            results.append(row)
            issues.append({"blob": e["git_blob_sha"], "issue": "unclassified unreadable TIFF; completeness fails", "error": str(exc)})
            continue
        with dataset as s:
            row["shape"] = list(s.shape); row["bands"] = s.count
            row["exact_grid"] = s.shape == shape and s.crs == crs and s.transform == transform and s.count == 1
            if s.count != 1:
                issues.append({"blob": e["git_blob_sha"], "issue": "unexpected multi-band prior output"}); continue
            raw = s.read(1)
            # A legacy nodata=0 tag must not hide almost the entire *raw*
            # binary prediction. Zero is a legitimate confidence value. Keep
            # those stored zeros for comparison rather than selecting only 1s.
            real = np.isfinite(raw) & ((s.read_masks(1) > 0) | ((raw == 0) & (s.nodata == 0)))
            row["legacy_nodata_zero_compared_as_raw_zero"] = bool(s.nodata == 0)
            if s.shape == shape:
                mask = foot & real
                row["comparison_mode"] = "raw pixel-index correspondence; no clipping or probability transform"
                row["comparison_coverage"] = float(mask.sum() / foot.sum())
                row["invalid_in_template_footprint"] = int(np.count_nonzero(foot & ~real))
                row["canonical_sha256"] = canonical_pixel_sha256(raw, ~foot)
                row["canonical_hash_equal"] = row["canonical_sha256"] == canonical
                if mask.any():
                    row.update(compare_vectors(candidate[mask], raw[mask]))
                else:
                    issues.append({"blob": e["git_blob_sha"], "issue": "no valid common pixels"}); continue
                if row["comparison_coverage"] < .90:
                    # A partial view alone is insufficient. Only these exact
                    # byte-pinned variants can be bridged, and only by proving
                    # their relationship to an actually published full field.
                    counterpart = NAN_STORAGE_COUNTERPARTS.get(raw_sha)
                    proved = False
                    if counterpart and any(item["git_blob_sha"] == counterpart[0] and item["sha256"] == counterpart[1] for item in entries):
                        companion = cache / f"{counterpart[0]}.tif"
                        if companion.is_file() and file_sha256(companion) == counterpart[1]:
                            with rasterio.open(companion) as complete_ds:
                                same = (complete_ds.count == 1 and complete_ds.shape == s.shape
                                        and complete_ds.crs == s.crs and complete_ds.transform == s.transform)
                                other = complete_ds.read(1)
                            proved = same and nan_zero_storage_equivalent(raw, other)
                            if proved:
                                row["verified_storage_counterpart"] = {
                                    "git_blob_sha": counterpart[0], "sha256": counterpart[1],
                                    "entire_array_equal_with_nan_represented_as_zero": True,
                                    "inside_missing_cells": int(np.count_nonzero(foot & np.isnan(raw))),
                                    "counterpart_nonzero_at_missing_cells": int(np.count_nonzero(other[foot & np.isnan(raw)])),
                                    "interpretation": "This row's Pearson uses ONLY its actual finite values. The separately published, byte-verified all-finite companion supplies the full-footprint comparison. No prior values were silently imputed."}
                    if not proved:
                        issues.append({"blob": e["git_blob_sha"], "issue": "less than 90% footprint raw comparison without an exact published storage counterpart"})
            else:
                row["canonical_hash_equal"] = False
                row["near_duplicate_h4"] = False
                row["comparison_mode"] = "different-size published format fixture; supplementary georeferenced comparison only"
            if not row["exact_grid"]:
                if s.crs is None:
                    issues.append({"blob": e["git_blob_sha"], "issue": "cannot align ungeoreferenced prior"}); continue
                aligned = np.full(shape, np.nan, dtype=np.float32)
                src = np.asarray(raw, dtype=np.float32); src[~real] = np.nan
                reproject(src, aligned, src_transform=s.transform, src_crs=s.crs,
                          src_nodata=np.nan, dst_transform=transform, dst_crs=crs,
                          dst_nodata=np.nan, resampling=Resampling.nearest, num_threads=2)
                common = foot & np.isfinite(aligned)
                coverage = float(common.sum() / foot.sum())
                supplemental = {"coverage": coverage, "interpretation": "partial spatial overlap is diagnostic, not a whole-map duplicate claim"}
                if common.any():
                    supplemental.update(compare_vectors(candidate[common], aligned[common]))
                    if coverage >= .9:
                        row["near_duplicate_h4"] |= supplemental["near_duplicate_h4"]
                row["aligned_comparison"] = supplemental
                if s.shape != shape and raw_sha != SPATIAL_FORMAT_FIXTURE_SHA256:
                    issues.append({"blob": e["git_blob_sha"], "issue": "different-size output is not the exact-byte known non-submission fixture"})
            row["near_duplicate_h4"] |= row["file_hash_equal"] or row["canonical_hash_equal"]
        results.append(row)
        if i % 25 == 0 or i == len(entries):
            print(f"full-resolution raw audit {i}/{len(entries)}", flush=True)
    near = [r for r in results if r.get("near_duplicate_h4")]
    full = [r for r in results if "pearson_r" in r]
    return {
        "candidate_sha256": sha, "candidate_canonical_sha256": canonical,
        "inventory_sha256": file_sha256(inventory_path), "rule": RULE,
        "scope": inventory.get("scope"), "inventory_blobs": len(entries),
        "hashed_and_audited_blobs": len(results), "full_size_raw_comparisons": len(full),
        "nonprediction_artifacts_forensically_classified": sum(r.get("nonprediction_artifact", False) for r in results),
        "low_coverage_storage_variants_with_verified_full_counterparts": sum("verified_storage_counterpart" in r for r in results),
        "raw_comparisons_covering_at_least_90pct_footprint": sum(r.get("comparison_coverage", 0) >= .9 for r in results),
        "nonprediction_policy": "The exact 110-byte header-only malformed artifact and 138-byte one-pixel demo are hash-accounted, not falsely correlated. Every other unreadable/unclassified artifact still fails completeness.",
        "exact_grid_comparisons": sum(r["exact_grid"] for r in results),
        "supplementary_aligned_comparisons": sum("aligned_comparison" in r for r in results),
        "completeness_pass": bool(full) and not issues and len(results) == len(entries),
        "uniqueness_pass": bool(full) and not near and not issues and len(results) == len(entries),
        "near_duplicate_count": len(near),
        "historic_h2b_rule_flags": sum(r.get("historic_h2b_rule_near_duplicate", False) for r in results),
        "max_abs_pearson": max((abs(r["pearson_r"]) for r in full), default=None),
        "max_top_jaccard": max((r["top_jaccard"] for r in full), default=None),
        "max_top_containment": max((r["top_containment"] for r in full), default=None),
        "closest_by_pearson": sorted(full, key=lambda r: abs(r["pearson_r"]), reverse=True)[:5],
        "issues": issues, "comparisons": results,
    }
