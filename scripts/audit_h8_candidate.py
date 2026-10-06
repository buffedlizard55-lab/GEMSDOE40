#!/usr/bin/env python3
"""Independent audit of the H8 candidate: format, novelty, and named surrogates.

Reads the generated candidate from ignored ``work/`` and the cached prior corpus.
Publishes nothing unless every format gate passes and no prior raster is a
near-duplicate.  Every surrogate number is labelled with the truth it used; none
of them is an organizer score.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np
import rasterio
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gemsdoe40.emission import credit_of_binary, pred_credit_of_binary  # noqa: E402
from gemsdoe40.raster import canonical_pixel_sha256, validate_candidate  # noqa: E402

BLOCKS = (4, 6)
GUARD = 3
PUBLISH_MAX_ABS_PEARSON = 0.95
PUBLISH_MAX_JACCARD = 0.50


def load(path: Path, band: int = 1) -> np.ndarray:
    with rasterio.open(path) as src:
        return src.read(band)


def credit_fast(prediction: np.ndarray, truth: np.ndarray) -> tuple[float, float, float]:
    """(TPw, C_pred, mass) via distance transforms -- O(N) instead of 49 shifts.

    math identical to gemsdoe40.emission.credit_of_binary / pred_credit_of_binary,
    which the test suite pins against the published equations.
    """
    p = np.asarray(prediction, dtype=np.float64)
    g = np.asarray(truth, dtype=bool)
    mass = float(p.sum())
    if not g.any():
        return 0.0, 0.0, mass
    d_g = ndimage.distance_transform_edt(~g)
    k_g = np.maximum(1.0 - d_g / 3.0, 0.0)
    # TPw = sum over TRUTH pixels of the best kernel credit nearby
    #     = sum over truth pixels of max over the 5x5..7x7 neighbourhood of p*k
    best = np.zeros_like(p)
    for dy in range(-3, 4):
        for dx in range(-3, 4):
            d = float(np.hypot(dx, dy))
            if d > 3.0:
                continue
            shifted = np.zeros_like(p)
            ys = slice(max(0, dy), g.shape[0] + min(0, dy))
            yd = slice(max(0, -dy), g.shape[0] + min(0, -dy))
            xs = slice(max(0, dx), g.shape[1] + min(0, dx))
            xd = slice(max(0, -dx), g.shape[1] + min(0, -dx))
            shifted[ys, xs] = p[yd, xd]
            np.maximum(best, shifted * (1.0 - d / 3.0), out=best)
    tp = float(best[g].sum())
    pred_credit = float((p * k_g).sum())
    return tp, pred_credit, mass


def stats_vs(truth: np.ndarray, prediction: np.ndarray, mask: np.ndarray | None = None) -> dict:
    pred = np.asarray(prediction, dtype=np.float32)
    g = np.asarray(truth, dtype=bool)
    if pred.shape != g.shape:
        raise ValueError("shape mismatch")
    if mask is not None:
        m = np.asarray(mask, dtype=bool)
        pred = np.where(m, 0.0, pred)
        g = g & ~m
    tp, pc, mass = credit_fast(pred, g)
    n = float(g.sum())
    fn = n - tp
    fp = mass - pc
    denom = tp + 0.2 * fp + 0.8 * fn
    return {"mass": mass, "tp": tp, "pred_credit": pc, "fp": fp, "fn": fn, "n_truth": n,
            "w": tp / mass if mass else 0.0, "cover": tp / n if n else 0.0,
            "dti_surrogate": tp / denom if denom > 0 else 0.0}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, default=ROOT / "work" / "h8" / "h8-euler-depthcluster-candidate.tif")
    parser.add_argument("--prior-cache", type=Path, default=ROOT / "data" / "prior")
    parser.add_argument("--publish", action="store_true", help="copy into docs/downloads when all gates pass")
    args = parser.parse_args()

    template = ROOT / "data" / "sample_submission.tif"
    receipt = validate_candidate(args.candidate, template)
    print(f"format valid: {receipt['valid']}  issues: {receipt['issues']}")
    if not receipt["valid"]:
        print("REFUSING: format gate failed")
        return 1

    with rasterio.open(template) as ds:
        footprint = np.isfinite(ds.read(1))
    pred = load(args.candidate)
    inside = pred[footprint]
    positive = inside > 0

    # ---- surrogates ------------------------------------------------------- #
    labels = load(ROOT / "data" / "labels.tif") == 1
    sgmc = load(ROOT / "data" / "external" / "derived_sgmc_faults_100m_u8.tif") == 1
    d_lab = ndimage.distance_transform_edt(~labels)
    sgmc_off = sgmc & (d_lab > 3)
    catalogue = {"truth": "data/labels.tif USGS+INGENIOUS catalogue (organizer-provided)",
                 "mask": "catalogue pixels excluded from the emitted domain (staff ruling, forum 11516/2)",
                 **stats_vs(labels, pred, mask=labels)}
    offcat = {"truth": "SGMC faults >300 m from the catalogue (owner-derived mirror, NOT organizer truth)",
              "mask": None, **stats_vs(sgmc_off, pred)}

    # ---- spatially blocked (4x6 with 3-cell guards) ------------------------ #
    rows, cols = footprint.shape
    bh, bw = rows // BLOCKS[0], cols // BLOCKS[1]
    blocks = []
    for i in range(BLOCKS[0]):
        for j in range(BLOCKS[1]):
            r0, r1 = i * bh, (i + 1) * bh if i < BLOCKS[0] - 1 else rows
            c0, c1 = j * bw, (j + 1) * bw if j < BLOCKS[1] - 1 else cols
            rs, re = min(r0 + GUARD, r1), max(r1 - GUARD, r0)
            cs, ce = min(c0 + GUARD, c1), max(c1 - GUARD, c0)
            blk = np.zeros_like(pred, dtype=bool)
            blk[rs:re, cs:ce] = True
            blk &= footprint
            n_cat = int((labels & blk).sum())
            n_off = int((sgmc_off & blk).sum())
            r0c, r1c = max(0, rs - 3), min(rows, re + 3)
            c0c, c1c = max(0, cs - 3), min(cols, ce + 3)
            crop = (slice(r0c, r1c), slice(c0c, c1c))
            entry = {"block": [i, j], "cells": int(blk.sum()), "catalogue_cells": n_cat,
                     "sgmc_offcat_cells": n_off, "emitted_cells": int((pred > 0)[blk].sum())}
            if n_cat:
                entry["catalogue"] = stats_vs(labels[crop], pred[crop], mask=labels[crop])
            if n_off:
                entry["sgmc_offcat"] = stats_vs(sgmc_off[crop], pred[crop])
            blocks.append(entry)
    truth_blocks = [b for b in blocks if b["catalogue_cells"] > 0]
    off_blocks = [b for b in blocks if b["sgmc_offcat_cells"] > 0]

    # ---- novelty against the cached prior corpus --------------------------- #
    priors = []
    pred_flat = pred.astype(np.float64).ravel()
    pred_bin = (pred > 0)
    topk = int(positive.sum())
    order = np.argsort(-pred_flat, kind="stable")[:topk]
    top_mask = np.zeros(pred_flat.size, dtype=bool)
    top_mask[order] = True
    for path in sorted(args.prior_cache.glob("*.tif")):
        try:
            with rasterio.open(path) as src:
                if src.shape != pred.shape:
                    priors.append({"file": path.name, "comparable": False,
                                   "reason": f"shape {src.shape}"})
                    continue
                other = src.read(1).astype(np.float32)
        except Exception as exc:  # noqa: BLE001
            priors.append({"file": path.name, "comparable": False, "reason": f"unreadable: {exc}"[:80]})
            continue
        other = np.where(np.isfinite(other), other, 0.0)
        finite = footprint & np.isfinite(other)
        a = pred[finite].astype(np.float64)
        b = other[finite].astype(np.float64)
        if a.std() == 0 or b.std() == 0:
            pearson = 0.0
        else:
            pearson = float(np.corrcoef(a, b)[0, 1])
        ob = (other > 0)
        ob_flat = ob.ravel()
        inter = int((top_mask & ob_flat).sum())
        union = int((top_mask | ob_flat).sum())
        jaccard = inter / union if union else 0.0
        containment = inter / max(topk, 1)
        exact = bool(np.array_equal(pred, other, equal_nan=True))
        priors.append({"file": path.name, "comparable": True, "pearson": pearson,
                       "jaccard_topk": jaccard, "containment_topk": containment,
                       "other_positive": int(ob_flat.sum()), "exact_equal": exact})
    comparable = [p for p in priors if p.get("comparable")]
    max_abs = max((abs(p["pearson"]) for p in comparable), default=0.0)
    max_jac = max((p["jaccard_topk"] for p in comparable), default=0.0)
    max_cont = max((p["containment_topk"] for p in comparable), default=0.0)
    top_matches = sorted(comparable, key=lambda p: -abs(p["pearson"]))[:5]
    novel = (max_abs <= PUBLISH_MAX_ABS_PEARSON and max_jac <= PUBLISH_MAX_JACCARD
             and not any(p["exact_equal"] for p in comparable))

    document = {
        "candidate": receipt,
        "emitted_pixels": int(positive.sum()),
        "distinct_values": int(np.unique(inside[positive]).size) if positive.any() else 0,
        "value_min": float(inside[positive].min()) if positive.any() else None,
        "value_max": float(inside[positive].max()) if positive.any() else None,
        "catalogue_surrogate": catalogue,
        "offcat_sgmc_surrogate": offcat,
        "blocks": blocks,
        "truth_bearing_blocks": len(truth_blocks),
        "offcat_truth_blocks": len(off_blocks),
        "novelty": {
            "corpus_size": len(comparable),
            "unreadable_or_incomparable": len(priors) - len(comparable),
            "max_abs_pearson": max_abs, "max_jaccard_topk": max_jac,
            "max_containment_topk": max_cont,
            "thresholds": {"max_abs_pearson": PUBLISH_MAX_ABS_PEARSON,
                           "max_jaccard_topk": PUBLISH_MAX_JACCARD},
            "top_matches": top_matches, "novel": bool(novel),
        },
        "is_a_score": False,
        "caveats": [
            "The catalogue surrogate cannot rank a candidate that deliberately avoids the catalogue: "
            "the family's best-scoring published file (37,654 dots, owner-reported 0.2778) has catalogue w = 0.0074.",
            "The SGMC-off-catalogue surrogate is owner-derived, not organizer truth, and the 2026-10-06 "
            "instrument audit measured no ranking power for it (Spearman +0.68 overall, sign inversion at matched mass).",
            "No hidden label, no organizer score and no authenticated DrivenData session is used here.",
        ],
    }
    out = ROOT / "docs" / "data" / "h8-audit.json"
    out.write_text(json.dumps(document, indent=2) + "\n")
    print(f"emitted {document['emitted_pixels']} px; values {document['value_min']}..{document['value_max']}")
    print(f"catalogue surrogate: w={catalogue['w']:.4f} cover={catalogue['cover']:.4f} dti={catalogue['dti_surrogate']:.5f}")
    print(f"offcat-sgmc surrogate: w={offcat['w']:.4f} cover={offcat['cover']:.4f} dti={offcat['dti_surrogate']:.5f}")
    print(f"novelty: max|pearson|={max_abs:.4f} max jaccard={max_jac:.4f} novel={novel}")

    if args.publish:
        if not novel:
            print("REFUSING to publish: near-duplicate of a prior raw output")
            return 1
        digest12 = receipt["canonical_pixels_sha256"][:12]
        name = f"gemsdoe40-h8-euler-lineament-depthcluster-20261006-{digest12}.tif"
        destination = ROOT / "docs" / "downloads" / name
        shutil.copyfile(args.candidate, destination)
        document["published"] = {"path": str(destination.relative_to(ROOT)), "bytes": destination.stat().st_size,
                                 "sha256": hashlib.sha256(destination.read_bytes()).hexdigest()}
        out.write_text(json.dumps(document, indent=2) + "\n")
        print(f"published {destination}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
