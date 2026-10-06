"""Metric-optimal emission for the DOE GEMS distance-weighted Tversky index.

Exact identity (derived from the organizer's definitions, page 967):

    DTI = TPw / (TPw + alpha*FPw + beta*FNw),   FNw = |G| - TPw
        = TPw / (alpha*(TPw + FPw) + beta*|G|)
        = TPw / (0.2*(TPw + FPw) + 0.8*|G|)          (alpha=0.2, beta=0.8)

Two emission consequences that this module implements:

1.  *Sparse* mass along a trace is far more mass-efficient than contiguous
    mass.  One dot on a 1-px trace is the arg-max for every truth cell within
    R=3 px, so it collects sum_d k(d) ~ 2.5-3.0 units of TPw for 1 unit of
    mass; a filled trace collects 1.0 unit of TPw per unit of mass.  Sparse
    mass therefore buys the same TPw for roughly a third of the FP budget when
    the trace hypothesis is right, and costs a third of the FP budget when it
    is wrong (measured in docs/data/emission_theory.json).
2.  The marginal rule.  Adding one unit of mass with expected TPw gain ``t``
    and expected FPw cost ``d`` changes DTI by

        gain  <=>  t * (0.2*FPw + 0.8*|G|)  >  0.2 * TPw * (t + d)

    which, for a dot that is either on the trace or useless (d = 1 - t),
    reduces to a very low bar: with DTI ~ 0.3 a dot needs only a few percent
    chance of landing on a hidden trace to be worth emitting.  The emitter
    therefore should *not* be pruned aggressively; it should be pruned exactly
    where the hidden-fault probability is ~0 (and never on top of the known
    catalogue, which the organizers exclude from scoring).
"""
from __future__ import annotations

import numpy as np


def nms_dots(
    field: np.ndarray,
    footprint: np.ndarray,
    *,
    min_separation_px: int = 3,
    budget: int | None = None,
    catalogue: np.ndarray | None = None,
    catalogue_buffer_px: int = 2,
    tie_break_row_major: bool = True,
) -> np.ndarray:
    """Greedy non-maximum-suppressed dot emission with a hard minimum separation.

    Cells are visited in descending field value (row-major for exact ties) and
    accepted only when no already-accepted dot lies within
    ``min_separation_px`` (Chebyshev distance, so 3 px means a 7x7 exclusion
    box).  The emitted raster is binary 1.0 at accepted dots and 0.0 elsewhere;
    outside ``footprint`` it is 0.0 (the writer sets NaN).
    """
    if min_separation_px < 1:
        raise ValueError("min_separation_px must be >= 1")
    h, w = field.shape
    positive = np.flatnonzero(field.ravel() > 0.0)
    if positive.size == 0:
        return np.zeros(field.shape, dtype=np.float32)
    # Stable descending sort keeps row-major order among exactly equal values.
    local = np.argsort(-field.ravel()[positive], kind="stable")
    order = positive[local]
    del tie_break_row_major  # documented behaviour: stable, row-major for ties
    blocked = np.zeros((h, w), dtype=bool)
    if catalogue is not None and catalogue_buffer_px >= 0:
        from scipy.ndimage import binary_dilation

        cat = np.asarray(catalogue, dtype=bool) & footprint
        if catalogue_buffer_px > 0 and cat.any():
            cat = binary_dilation(cat, iterations=int(catalogue_buffer_px))
        blocked |= cat
    blocked |= ~np.asarray(footprint, dtype=bool)

    out = np.zeros((h, w), dtype=np.float32)
    sep = int(min_separation_px)
    # Exclusion box half-width: a dot closer than `sep` is too close, so the
    # exclusion radius is sep-1 in Chebyshev terms.
    excl = max(sep - 1, 0)
    accepted = 0
    candidate_mask = np.zeros((h, w), dtype=bool)
    flat_blocked = blocked.ravel()
    for idx in order:
        if flat_blocked[idx]:
            continue
        r, c = divmod(int(idx), w)
        if candidate_mask[r, c]:
            continue
        out[r, c] = 1.0
        accepted += 1
        r0, r1 = max(0, r - excl), min(h, r + excl + 1)
        c0, c1 = max(0, c - excl), min(w, c + excl + 1)
        candidate_mask[r0:r1, c0:c1] = True
        for rr in range(r0, r1):
            flat_blocked[rr * w + c0: rr * w + c1] = True
        if budget is not None and accepted >= int(budget):
            break
    return out


def precision_curve(
    field: np.ndarray,
    truth: np.ndarray,
    footprint: np.ndarray,
    *,
    budgets: tuple[int, ...],
    min_separation_px: int = 3,
    catalogue: np.ndarray | None = None,
    catalogue_buffer_px: int = 2,
) -> list[dict]:
    """Emitted mass, TPw, FPw and DTI for a ladder of dot budgets.

    Computed on the *whole* grid; callers are responsible for holding out
    spatial blocks when they use this ladder to choose a budget.
    """
    from .research_metric import metric_components

    rows = []
    for b in budgets:
        dots = nms_dots(
            field, footprint, min_separation_px=min_separation_px, budget=int(b),
            catalogue=catalogue, catalogue_buffer_px=catalogue_buffer_px,
        )
        comp = metric_components(dots, truth, valid=footprint)
        rows.append({
            "budget": int(b),
            "emitted_px": int((dots > 0).sum()),
            "tpw": comp.tp,
            "fpw": comp.fp,
            "fnw": comp.fn,
            "truth_px": comp.n_truth,
            "dti": comp.score,
            "credit_per_dot": float(comp.tp / max(int((dots > 0).sum()), 1)),
        })
    return rows
