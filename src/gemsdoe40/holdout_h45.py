#!/usr/bin/env python3
"""A non-circular, spatially-blocked validation instrument for fault candidates.

Why the inherited instrument cannot be used
-------------------------------------------
The historical instrument in this repository scored candidates against an SGMC-derived
fault raster.  Two defects were already recorded there and are confirmed again here:

1. *Circularity*: the incumbent comparator is itself derived from that same catalogue, and
   the "truth" is the catalogue, so agreement with the catalogue is rewarded even though
   the organiser **masks all known faults out of scoring**
   (https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516).
2. *Inverted relationship*: across 32 prior predictions with owner-reported public scores,
   the Spearman correlation between catalogue hits and the reported score is **negative**
   (rho = -0.66 for hits within 300 m of ``labels.tif``, -0.79 against the SGMC raster).
   Predicting the catalogue is therefore not merely uninformative, it is harmful.

What this instrument does instead
---------------------------------
For each block of a K x L partition:

* the mapped faults **inside** the block are treated as if they had never been mapped
  (this is the entire premise of the competition: the test labels are faults that are
  absent from the public catalogue);
* a guard band around the block is deleted from both truth and prediction, so a fault
  crossing the block boundary is not scored as a false positive just outside the block;
* the official distance-weighted Tversky index is computed with alpha = 0.2, beta = 0.8
  and the 300 m triangular kernel, over the whole footprint, with that block's faults as
  the only truth.

The score is then averaged over blocks.  Because the candidate is built from geophysics
only and never reads the label raster, treating one block's faults as hidden is a genuine
out-of-sample test of "would this method have found an unmapped fault here?".

Two companion diagnostics are reported:

* ``lift_blocked``: the top-``budget`` pixel hit rate against mapped faults, divided by the
  footprint base rate, averaged over blocks.  A value of 1.0 means no skill.
* ``catalogue_gradient``: hit *density* of the candidate as a function of distance to the
  mapped catalogue, pooled over blocks.  The owner's own submissions show that pixels
  within 200 m of the catalogue carry far less real credit than the average
  (0.2708 -> 0.2778 after deleting 2,545 such pixels out of 40,199), so a candidate that
  dumps mass into the catalogue halo is measurably worse.

Neither statistic uses hidden labels, an organiser receipt, or a prior submission's pixels.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Iterable, Sequence

import numpy as np
from scipy.ndimage import binary_dilation, binary_erosion, distance_transform_edt

ALPHA = 0.2
BETA = 0.8
RADIUS_PX = 3.0
EPS = 1e-12

_OFFSETS: tuple[tuple[int, int, float], ...] = tuple(
    (dy, dx, max(1.0 - math.hypot(dy, dx) / RADIUS_PX, 0.0))
    for dy in range(-3, 4) for dx in range(-3, 4)
    if max(1.0 - math.hypot(dy, dx) / RADIUS_PX, 0.0) > 0.0
)


def triangular_credit(pred: np.ndarray, truth: np.ndarray,
                      active: np.ndarray | None = None) -> dict:
    """Official distance-weighted TP_w, FP_w, FN_w for a soft prediction raster.

    Identical algebra to ``src/gems40/metric.py``; duplicated (not imported) so that this
    instrument can be validated independently against the organiser's published equations.
    """
    active = np.ones(pred.shape, bool) if active is None else active
    p = np.where(active & np.isfinite(pred), np.asarray(pred, dtype=np.float64), 0.0)
    p = np.maximum(p, 0.0)
    g = active & (np.asarray(truth) > 0)
    n_g = int(g.sum())
    if n_g == 0:
        return dict(dti=0.0, tp_w=0.0, fp_w=float(p.sum()), fn_w=0.0, n_g=0)
    yy, xx = np.nonzero(g)
    credit = np.zeros(yy.size, dtype=np.float64)
    for dy, dx, k in _OFFSETS:
        ny, nx = yy + dy, xx + dx
        ok = (ny >= 0) & (ny < pred.shape[0]) & (nx >= 0) & (nx < pred.shape[1])
        credit[ok] = np.maximum(credit[ok], p[ny[ok], nx[ok]] * k)
    tp_w = float(credit.sum())
    near = np.maximum(1.0 - distance_transform_edt(~g) / RADIUS_PX, 0.0)
    fp_w = float((p * (1.0 - near)).sum())
    fn_w = float(n_g) - tp_w
    dti = tp_w / (tp_w + ALPHA * fp_w + BETA * fn_w + EPS)
    return dict(dti=float(dti), tp_w=tp_w, fp_w=fp_w, fn_w=fn_w, n_g=n_g,
                mass=float(p.sum()), support=int((p > 0).sum()))


@dataclass
class BlockedHoldout:
    rows: int = 6
    cols: int = 5
    guard_px: int = 12
    budget_px: int = 40000
    blocks: list[tuple[int, int, int, int]] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.blocks = []


def partition(shape: tuple[int, int], rows: int, cols: int) -> list[tuple[int, int, int, int]]:
    h, w = shape
    by, bx = h // rows, w // cols
    return [(r * by, (r + 1) * by if r < rows - 1 else h,
             c * bx, (c + 1) * bx if c < cols - 1 else w)
            for r in range(rows) for c in range(cols)]


def evaluate(candidate: np.ndarray, truth_full: np.ndarray, valid: np.ndarray,
             rows: int = 6, cols: int = 5, guard_px: int = 12,
             budget_px: int = 40000) -> dict:
    """Blocked holdout DTI plus the two companion diagnostics."""
    pred = np.where(valid, np.asarray(candidate, dtype=np.float64), 0.0)
    pred = np.maximum(np.where(np.isfinite(pred), pred, 0.0), 0.0)
    blocks = partition(valid.shape, rows, cols)
    per_block = []
    for y0, y1, x0, x1 in blocks:
        inside = np.zeros(valid.shape, bool)
        inside[y0:y1, x0:x1] = True
        # guard: erode the block, then delete a band around it from truth and prediction
        if guard_px:
            kept = binary_erosion(inside & valid, np.ones((2 * guard_px + 1,) * 2, bool))
        else:
            kept = inside & valid
        truth_b = truth_full & kept
        if truth_b.sum() < 50:
            continue
        scored = valid & ~binary_dilation(inside & ~kept, np.ones((2 * guard_px + 1,) * 2, bool))
        res = triangular_credit(pred, truth_b, scored)
        res["block"] = (int(y0), int(y1), int(x0), int(x1))
        res["truth_px"] = int(truth_b.sum())
        per_block.append(res)

    dti_mean = float(np.mean([b["dti"] for b in per_block])) if per_block else float("nan")
    dti_median = float(np.median([b["dti"] for b in per_block])) if per_block else float("nan")

    # ---- blocked lift -----------------------------------------------------------------
    d_cat = distance_transform_edt(~truth_full)
    halo = (d_cat <= RADIUS_PX) & valid
    base = halo.sum() / max(valid.sum(), 1)
    lifts = []
    for y0, y1, x0, x1 in blocks:
        sub_valid = valid[y0:y1, x0:x1]
        sub_halo = halo[y0:y1, x0:x1]
        if sub_valid.sum() < 2000 or sub_halo.sum() < 50:
            continue
        budget = max(50, int(budget_px * sub_valid.sum() / max(valid.sum(), 1)))
        sub = np.where(sub_valid, pred[y0:y1, x0:x1], -np.inf).ravel()
        idx = np.argpartition(sub, -budget)[-budget:]
        hit = sub_halo.ravel()[idx].mean()
        base_b = sub_halo.sum() / sub_valid.sum()
        lifts.append(float(hit / base_b) if base_b > 0 else float("nan"))

    # ---- catalogue-distance gradient ---------------------------------------------------
    flat_pred = pred.ravel()
    flat_d = d_cat.ravel()
    flat_valid = valid.ravel()
    bands = [(0, 1), (1, 2), (2, 3), (3, 5), (5, 10), (10, 1 << 30)]
    gradient = {}
    total_mass = float(flat_pred[flat_valid].sum())
    for lo, hi in bands:
        sel = flat_valid & (flat_d >= lo) & (flat_d < hi)
        area = int(sel.sum())
        mass = float(flat_pred[sel].sum())
        gradient[f"{lo*100:g}-{hi*100:g}m"] = dict(
            area=area, mass=mass,
            mass_fraction=float(mass / total_mass) if total_mass else 0.0,
            density=float(mass / area) if area else 0.0)
    overall_density = total_mass / max(int(flat_valid.sum()), 1)
    for k in gradient:
        gradient[k]["relative_density"] = (gradient[k]["density"] / overall_density
                                           if overall_density > 0 else 0.0)

    return dict(
        dti_blocked_mean=dti_mean, dti_blocked_median=dti_median,
        blocks_used=len(per_block), per_block=per_block,
        lift_blocked_mean=float(np.mean(lifts)) if lifts else float("nan"),
        lift_blocked_median=float(np.median(lifts)) if lifts else float("nan"),
        lift_blocks_used=len(lifts),
        catalogue_gradient=gradient,
        total_mass=total_mass,
        support=int((pred > 0).sum()),
    )


# --------------------------------------------------------------------------------------
# GAP-BLOCK holdout
# --------------------------------------------------------------------------------------
# Rationale (measured, not assumed): across the owner's own submissions, agreement with the
# public catalogue is monotonically ANTI-correlated with the reported public score - the
# 0.2600 submission has a catalogue lift of 2.27 while the 0.2778 submission has 0.70,
# because the latter deletes every prediction within 200 m of the catalogue.  The natural
# reading is that the organiser's expert-labelled test faults lie in areas the public
# catalogue does not already cover, so the useful question is not "does this candidate find
# mapped faults?" but "does it find faults IN CATALOGUE GAPS?".
#
# Construction, per block (block masks are rectangles, so all set algebra is O(1)):
#   * the block's own faults are treated as unmapped (the competition premise);
#   * a guard band around the block is removed from both truth and the scored region;
#   * truth is further restricted to pixels at least gap_px from every *remaining*
#     (still-mapped) fault, so the truth really does sit in a catalogue gap;
#   * the official DTI is computed over the whole footprint against that truth, which
#     penalises mass dumped onto still-mapped catalogue faults.


class GapHoldout:
    """Precomputes the per-block truth / gap / scored masks once, then scores candidates."""

    def __init__(self, truth_full: np.ndarray, valid: np.ndarray, rows: int = 5,
                 cols: int = 4, guard_px: int = 12, gap_px: int = 10,
                 budget_px: int = 40000, min_truth_px: int = 50):
        self.valid = valid
        self.budget_px = budget_px
        self.blocks = []
        h, w = valid.shape
        by, bx = h // rows, w // cols
        for r in range(rows):
            for c in range(cols):
                y0, y1 = r * by, (r + 1) * by if r < rows - 1 else h
                x0, x1 = c * bx, (c + 1) * bx if c < cols - 1 else w
                g = guard_px
                ky0, ky1 = max(y0 + g, 0), min(y1 - g, h)          # truth (eroded) window
                kx0, kx1 = max(x0 + g, 0), min(x1 - g, w)
                py0, py1 = max(y0 - g, 0), min(y1 + g, h)          # padded block
                px0, px1 = max(x0 - g, 0), min(x1 + g, w)
                sy0, sy1 = max(y0 + 2 * g, 0), min(y1 - 2 * g, h)  # scored interior strip
                sx0, sx1 = max(x0 + 2 * g, 0), min(x1 - 2 * g, w)
                known_rest = truth_full.copy()
                known_rest[py0:py1, px0:px1] = False
                dist_known = distance_transform_edt(~known_rest)
                gap = dist_known > gap_px
                truth_b = np.zeros(valid.shape, bool)
                truth_b[ky0:ky1, kx0:kx1] = truth_full[ky0:ky1, kx0:kx1] & gap[ky0:ky1, kx0:kx1]
                scored = valid.copy()
                scored[py0:py1, px0:px1] = False
                scored[sy0:sy1, sx0:sx1] = True                     # reinstate block interior
                scored &= valid
                # gap-lift window: the block, restricted to the gap
                sub_valid = np.zeros(valid.shape, bool)
                sub_valid[y0:y1, x0:x1] = True
                sub_valid &= valid & gap
                sub_truth = np.zeros(valid.shape, bool)
                sub_truth[y0:y1, x0:x1] = truth_full[y0:y1, x0:x1]
                sub_truth &= gap
                if int(truth_b.sum()) < min_truth_px:
                    continue
                halo = (distance_transform_edt(~sub_truth) <= RADIUS_PX) & sub_valid
                area = int(sub_valid.sum())
                base = float(halo.sum() / area) if area else 0.0
                budget = max(50, int(budget_px * area / max(int(valid.sum()), 1)))
                self.blocks.append(dict(rect=(y0, y1, x0, x1), truth=truth_b,
                                        scored=scored, sub_valid=sub_valid,
                                        halo=halo, base=base, budget=budget,
                                        truth_px=int(truth_b.sum())))
        self.total_truth_px = int(sum(b["truth_px"] for b in self.blocks))

    def __call__(self, candidate: np.ndarray) -> dict:
        pred = np.where(self.valid, np.asarray(candidate, dtype=np.float64), 0.0)
        pred = np.maximum(np.where(np.isfinite(pred), pred, 0.0), 0.0)
        dtis, lifts, pb = [], [], []
        for b in self.blocks:
            res = triangular_credit(pred, b["truth"], b["scored"])
            dtis.append(res["dti"])
            pb.append(dict(dti=res["dti"], tp_w=res["tp_w"], fp_w=res["fp_w"],
                           fn_w=res["fn_w"], truth_px=b["truth_px"], rect=b["rect"]))
            (y0, y1, x0, x1) = b["rect"]
            sub = np.where(b["sub_valid"][y0:y1, x0:x1], pred[y0:y1, x0:x1], -np.inf).ravel()
            idx = np.argpartition(sub, -b["budget"])[-b["budget"]:]
            hit = float(b["halo"][y0:y1, x0:x1].ravel()[idx].mean())
            lifts.append(hit / b["base"] if b["base"] > 0 else float("nan"))
        return dict(
            gap_dti_mean=float(np.mean(dtis)) if dtis else float("nan"),
            gap_dti_median=float(np.median(dtis)) if dtis else float("nan"),
            gap_blocks_used=len(self.blocks),
            gap_truth_px=self.total_truth_px,
            gap_lift_mean=float(np.nanmean(lifts)) if lifts else float("nan"),
            gap_lift_median=float(np.nanmedian(lifts)) if lifts else float("nan"),
            gap_lift_blocks=len(lifts),
            per_block=pb)
