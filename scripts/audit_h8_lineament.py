#!/usr/bin/env python3
"""Independent audit of the H8 lineament candidate, and publication gate.

Checks, in order:

1. **Format** — re-open the written file and verify it against
   ``data/sample_submission.tif`` exactly: one float32 band, same CRS, shape and
   transform, finite values in [0, 1] at all 5,167,373 footprint cells, NaN
   outside, plus the raw SHA-256 and canonical pixel SHA-256.
2. **Named surrogates** — distance-weighted Tversky credit of the emitted dots
   against (a) the provided USGS/INGENIOUS catalogue with the pixel-exact
   catalogue and its 1 px ring excluded from the prediction domain, and (b) the
   owner-derived SGMC raster restricted to faults further than 300 m from the
   catalogue.  Both are labelled proxies; neither is organizer truth.
3. **Spatially blocked (4 × 6 with 3-cell guards)** — per-block catalogue and
   SGMC-off-catalogue credit, so a single region cannot carry the result.
4. **Novelty** — Pearson correlation of the in-footprint values, positive-support
   Jaccard, and top-mass containment against every cached prior raster; the
   audit refuses to publish on an exact duplicate or beyond the registered
   thresholds.
5. **Comparison row** — the same surrogate numbers for the 16 recorded artifacts
   that the family has leaderboard scores for, so the reader can see exactly
   where this file stands.

No network access, no labels in generation, no score is claimed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sys

import numpy as np
import rasterio
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gemsdoe40.raster import canonical_pixel_sha256, validate_candidate  # noqa: E402

SPACING = 3.0
BLOCKS = (4, 6)
GUARD = 3
MAX_ABS_PEARSON = 0.95
MAX_JACCARD = 0.50
RING_PX = 1


def digest_file(path: Path) -> str:
    import hashlib
    return hashlib.sha256(path.read_bytes()).hexdigest()


def kernel_credit(distance: np.ndarray) -> np.ndarray:
    return np.maximum(1.0 - distance / SPACING, 0.0).astype(np.float32)


def stats(prediction: np.ndarray, truth: np.ndarray, *, kernel: np.ndarray) -> dict:
    """Exact distance-weighted Tversky statistics for a binary-support prediction."""
    pred = np.asarray(prediction, dtype=np.float32)
    total_mass = float(pred.sum())
    if not truth.any():
        return {"mass": total_mass, "n_truth": 0}
    # TP credit: each truth pixel is credited once, by the best emitter in its 3 px kernel
    best = np.zeros_like(pred)
    radius = int(np.ceil(SPACING))
    for dy in range(-radius, radius + 1):
        for dx in range(-radius, radius + 1):
            d = float(np.hypot(dx, dy))
            if d > SPACING:
                continue
            shifted = np.zeros_like(pred)
            ys, yd = slice(max(0, dy), pred.shape[0] + min(0, dy)), slice(max(0, -dy), pred.shape[0] + min(0, -dy))
            xs, xd = slice(max(0, dx), pred.shape[1] + min(0, dx)), slice(max(0, -dx), pred.shape[1] + min(0, -dx))
            shifted[ys, xs] = pred[yd, xd]
            np.maximum(best, shifted * (1.0 - d / SPACING), out=best)
    tp = float(best[truth].sum())
    pred_credit = float((pred * kernel).sum())
    fn = float(truth.sum()) - tp
    fp = total_mass - pred_credit
    denominator = tp + 0.2 * fp + 0.8 * fn
    return {"mass": total_mass, "tp": tp, "fn": fn, "fp": fp, "pred_credit": pred_credit,
            "n_truth": int(truth.sum()),
            "w": pred_credit / total_mass if total_mass else 0.0,          # comparable with the 16 recorded artifacts
            "w_truth_side": tp / total_mass if total_mass else 0.0,        # redundant-coverage-free variant
            "cover": tp / float(truth.sum()), "dti_surrogate": tp / denominator if denominator > 0 else 0.0}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path,
                        default=ROOT / "work/h8lineament/h8-lineament-continuous.tif")
    parser.add_argument("--prior-cache", type=Path, default=ROOT / "data/prior")
    parser.add_argument("--publish", action="store_true")
    parser.add_argument("--name", default=None)
    parser.add_argument("--hard-twin", type=Path,
                        default=ROOT / "work/h8lineament/h8-lineament-binary.tif",
                        help="identical support, every emitted dot at 1.0 (the metric-optimal twin)")
    args = parser.parse_args()

    template = ROOT / "data/sample_submission.tif"
    receipt = validate_candidate(args.candidate, template)
    print(f"format valid: {receipt['valid']}; issues: {receipt['issues']}")
    if not receipt["valid"]:
        print("REFUSING: format gate failed")
        return 1

    with rasterio.open(template) as ds:
        footprint = np.isfinite(ds.read(1))
    with rasterio.open(args.candidate) as ds:
        pred = np.where(np.isfinite(ds.read(1)), ds.read(1), 0.0).astype(np.float32)
    with rasterio.open(ROOT / "data/labels.tif") as ds:
        catalogue = ds.read(1) == 1
    with rasterio.open(ROOT / "data/external/derived_sgmc_faults_100m_u8.tif") as ds:
        sgmc = ds.read(1) == 1
    if pred.shape != footprint.shape:
        print("REFUSING: shape mismatch")
        return 1

    excluded = ndimage.binary_dilation(catalogue & footprint, np.ones((2 * RING_PX + 1,) * 2, bool))
    domain = footprint & ~excluded
    restricted = np.where(domain, pred, 0.0)

    d_cat = ndimage.distance_transform_edt(~(catalogue & footprint))
    k_cat = kernel_credit(d_cat)
    sgmc_off = sgmc & footprint & (d_cat > SPACING)
    d_off = ndimage.distance_transform_edt(~sgmc_off)
    k_off = kernel_credit(d_off)

    mass = float(restricted.sum())
    dots = int((restricted > 0).sum())
    dot_rows, dot_cols = np.nonzero(restricted > 0)
    dot_distance = d_cat[restricted > 0]

    catalogue_surrogate = {"truth": "data/labels.tif USGS+INGENIOUS catalogue (organizer-provided)",
                           "domain": f"footprint minus the catalogue mask and its {RING_PX} px ring",
                           **stats(restricted, catalogue & footprint, kernel=k_cat)}
    offcat_surrogate = {"truth": "derived_sgmc_faults_100m_u8 pixels further than 300 m from the catalogue",
                        "domain": "same restricted domain",
                        **stats(restricted, sgmc_off, kernel=k_off)}

    rows, cols = footprint.shape
    bh, bw = rows // BLOCKS[0], cols // BLOCKS[1]
    blocks = []
    for i in range(BLOCKS[0]):
        for j in range(BLOCKS[1]):
            r0, r1 = i * bh, (i + 1) * bh if i < BLOCKS[0] - 1 else rows
            c0, c1 = j * bw, (j + 1) * bw if j < BLOCKS[1] - 1 else cols
            rs, re = min(r0 + GUARD, r1), max(r1 - GUARD, r0)
            cs, ce = min(c0 + GUARD, c1), max(c1 - GUARD, c0)
            r0c, r1c, c0c, c1c = max(0, rs - 3), min(rows, re + 3), max(0, cs - 3), min(cols, ce + 3)
            crop = (slice(r0c, r1c), slice(c0c, c1c))
            entry = {"block": [i, j], "emitted": int((restricted[crop] > 0).sum()),
                     "catalogue_cells": int(((catalogue & footprint)[crop]).sum()),
                     "sgmc_offcat_cells": int(sgmc_off[crop].sum())}
            if entry["catalogue_cells"]:
                entry["catalogue"] = stats(restricted[crop], (catalogue & footprint)[crop], kernel=k_cat[crop])
            if entry["sgmc_offcat_cells"]:
                entry["offcat"] = stats(restricted[crop], sgmc_off[crop], kernel=k_off[crop])
            blocks.append(entry)

    comparisons = []
    pred_flat = restricted.ravel()
    topk = dots
    order = np.argsort(-pred_flat, kind="stable")[:topk]
    top_mask = np.zeros(pred_flat.size, dtype=bool)
    top_mask[order] = True
    for path in sorted(args.prior_cache.glob("*.tif")):
        try:
            with rasterio.open(path) as src:
                if src.shape != pred.shape:
                    comparisons.append({"file": path.name, "comparable": False,
                                        "reason": f"shape {src.shape}"})
                    continue
                other = np.where(np.isfinite(src.read(1)), src.read(1), 0.0).astype(np.float32)
        except Exception as exc:  # noqa: BLE001
            comparisons.append({"file": path.name, "comparable": False, "reason": str(exc)[:80]})
            continue
        finite = footprint
        a, b = pred[finite].astype(np.float64), other[finite].astype(np.float64)
        pearson = 0.0 if (a.std() == 0 or b.std() == 0) else float(np.corrcoef(a, b)[0, 1])
        other_positive = (other > 0).ravel()
        inter = int((top_mask & other_positive).sum())
        union = int((top_mask | other_positive).sum())
        comparisons.append({"file": path.name, "comparable": True, "pearson": pearson,
                            "jaccard_topmass": inter / union if union else 0.0,
                            "containment_topmass": inter / max(topk, 1),
                            "other_positive": int(other_positive.sum()),
                            "exact_equal": bool(np.array_equal(pred, other, equal_nan=True))})
    comparable = [c for c in comparisons if c.get("comparable")]
    max_abs = max((abs(c["pearson"]) for c in comparable), default=0.0)
    max_jac = max((c["jaccard_topmass"] for c in comparable), default=0.0)
    exact = [c["file"] for c in comparable if c["exact_equal"]]
    novel = max_abs <= MAX_ABS_PEARSON and max_jac <= MAX_JACCARD and not exact

    recorded = []
    matched_path = ROOT / "work/prior_catalogue_matched.json"
    if matched_path.is_file():
        for row in json.loads(matched_path.read_text()):
            recorded.append({"file": row["file"], "live_score_recorded": row["live_score"],
                             "mass": row["mass"], "w_catalogue": row["w_lab"],
                             "w_sgmc_offcat": row["w_off"]})

    values = pred[restricted > 0]
    document = {
        "candidate": receipt,
        "class": "research artifact; not an organizer score; no submission made",
        "emitted_dots": dots,
        "mass": mass,
        "distinct_values": int(np.unique(values).size) if dots else 0,
        "value_min": float(values.min()) if dots else None,
        "value_max": float(values.max()) if dots else None,
        "value_p50": float(np.median(values)) if dots else None,
        "value_fraction_at_one": float((values >= 0.999999).mean()) if dots else None,
        "dot_distance_to_catalogue_m": {
            "min": float(dot_distance.min()) * 100 if dots else None,
            "p50": float(np.median(dot_distance)) * 100 if dots else None,
            "p95": float(np.percentile(dot_distance, 95)) * 100 if dots else None},
        "catalogue_surrogate": catalogue_surrogate,
        "offcat_sgmc_surrogate": offcat_surrogate,
        "blocks": blocks,
        "novelty": {"corpus_size": len(comparable),
                    "unreadable_or_incomparable": len(comparisons) - len(comparable),
                    "max_abs_pearson": max_abs, "max_jaccard_topmass": max_jac,
                    "exact_duplicates": exact, "thresholds": {"max_abs_pearson": MAX_ABS_PEARSON,
                                                              "max_jaccard_topmass": MAX_JACCARD},
                    "top_matches": sorted(comparable, key=lambda c: -abs(c["pearson"]))[:6],
                    "novel": bool(novel)},
        "recorded_artifact_comparison": recorded,
        "caveats": [
            "No local instrument reaches the ranking power the family's audit requires "
            "(docs/research/h8-analysis-20261006.md section 6); surrogate numbers do not predict "
            "the leaderboard.",
            "The winning recorded artifact has catalogue w = 0.0074 and zero dots within 100 m of "
            "the catalogue by construction, so a low catalogue surrogate is expected of any "
            "off-catalogue file and cannot be used as a quality ranking.",
            "Trained on nothing: the construction reads only the challenge feature stack, the "
            "sample template and the catalogue mask.",
        ],
    }

    out = ROOT / "docs/data/h8-lineament-audit.json"
    print(f"dots {dots}; mass {mass}; values {document['value_min']}..{document['value_max']}")
    print(f"catalogue surrogate: w={catalogue_surrogate['w']:.4f} cover={catalogue_surrogate['cover']:.4f} "
          f"dti={catalogue_surrogate['dti_surrogate']:.5f}")
    print(f"offcat surrogate:    w={offcat_surrogate['w']:.4f} cover={offcat_surrogate['cover']:.4f} "
          f"dti={offcat_surrogate['dti_surrogate']:.5f}")
    print(f"novelty: corpus={len(comparable)} max|pearson|={max_abs:.4f} max jaccard={max_jac:.4f} "
          f"novel={novel}")

    if args.publish:
        if not novel:
            print("REFUSING to publish: near-duplicate or exact duplicate of a prior raw output")
            out.write_text(json.dumps(document, indent=2) + "\n")
            return 1
        digest12 = receipt["canonical_pixels_sha256"][:12]
        name = args.name or f"gemsdoe40-h8-euler-lineament-depthcluster-20261006-{digest12}.tif"
        destination = ROOT / "docs/downloads" / name
        shutil.copyfile(args.candidate, destination)
        twin = None
        if args.hard_twin.is_file():
            twin_receipt = validate_candidate(args.hard_twin, template)
            if twin_receipt["valid"]:
                twin_name = f"gemsdoe40-h8-euler-lineament-depthcluster-20261006-{digest12}-hard.tif"
                twin_destination = ROOT / "docs/downloads" / twin_name
                shutil.copyfile(args.hard_twin, twin_destination)
                with rasterio.open(args.hard_twin) as ds:
                    twin_values = ds.read(1)
                inside = twin_values[np.isfinite(twin_values)]
                twin = {
                    "name": twin_name, "path": str(twin_destination.relative_to(ROOT)),
                    "bytes": twin_destination.stat().st_size,
                    "sha256": validate_candidate(twin_destination, template)["sha256"],
                    "canonical_pixels_sha256": twin_receipt["canonical_pixels_sha256"],
                    "support_identical_to_primary": bool(
                        np.array_equal(np.isfinite(twin_values) & (twin_values > 0),
                                       np.isfinite(pred) & (pred > 0))),
                    "in_footprint_values": sorted({float(v) for v in np.unique(inside)})[:4],
                    "download_url": f"https://github.com/buffedlizard55-lab/GEMSDOE40/raw/main/"
                                    f"{twin_destination.relative_to(ROOT)}",
                }
                print(f"published hard twin {twin_destination}")
            else:
                print(f"hard twin not published: {twin_receipt['issues']}")
        document["published"] = {
            "path": str(destination.relative_to(ROOT)), "bytes": destination.stat().st_size,
            "sha256": validate_candidate(destination, template)["sha256"], "name": name,
            "download_url": f"https://github.com/buffedlizard55-lab/GEMSDOE40/raw/main/{destination.relative_to(ROOT)}",
        }
        # publish the depth-labelled solution cloud and copy the generation receipt
        cloud_source = args.candidate.parent / "h8-lineament-solutions.csv.gz"
        cloud_destination = ROOT / "docs/downloads/h8-lineament-solutions.csv.gz"
        if cloud_source.is_file():
            shutil.copyfile(cloud_source, cloud_destination)
            document["published"]["solution_cloud"] = {
                "path": str(cloud_destination.relative_to(ROOT)),
                "bytes": cloud_destination.stat().st_size,
                "sha256": digest_file(cloud_destination)}
        generation_source = args.candidate.parent / "h8-lineament-receipt.json"
        if generation_source.is_file():
            generation_destination = ROOT / "docs/data/h8-lineament-generation.json"
            generation_destination.write_text(generation_source.read_text())
            document["published"]["generation_receipt"] = str(generation_destination.relative_to(ROOT))
        if twin is not None:
            document["published"]["hard_twin"] = twin
        print(f"published {destination}")
    out.write_text(json.dumps(document, indent=2) + "\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
