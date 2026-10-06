"""Surrogate instruments that score a candidate the way the organiser scores it.

The organiser's published rule (DrivenData problem description, page 967, and the
staff clarification in community thread 11516) is:

* metric - distance-weighted Tversky index (DTI) with ``alpha = 0.2``,
  ``beta = 0.8`` and a triangular kernel of 300 m = 3 px support;
* ground truth - expert-identified faults that are **not** in the published
  USGS/INGENIOUS catalogue;
* the published catalogue itself is masked pixel-exactly, so prediction mass on
  a catalogue pixel is excluded from the false-positive term.

No instrument in this repository can know the hidden label set.  What can be
built from official inputs is a *surrogate* of the same shape: a truth set of
mapped structures that are absent from the published catalogue, plus the
pixel-exact catalogue mask.  Two surrogates are used here:

``sgmc_off_catalogue``
    USGS State Geologic Map Compilation faults (mrdata.usgs.gov) that are not
    catalogue pixels.  Large, but a different population: mostly older bedrock
    faults.  Beware circularity - a field built *from* SGMC is scored by
    construction; such fields are excluded from calibration.
``isolated_catalogue_components``
    Connected components of the published catalogue that sit at least
    ``min_separation_px`` from every other component.  Held out as "unseen
    faults" while the remaining catalogue is the mask.  Small, but
    non-gameable: the mask does not reveal where the held-out truth is, so a
    "paint the catalogue flank" strategy earns nothing here.

Both are scored with :func:`masked_dti`, which is the exact algebraic form of
the organiser's metric under the masking rule.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import binary_dilation, distance_transform_edt, label as cc_label

ALPHA = 0.2
BETA = 0.8
RADIUS_PX = 3.0
EPS = 1e-12


def kernel(d, radius: float = RADIUS_PX):
    return np.maximum(1.0 - np.asarray(d, dtype=np.float64) / radius, 0.0)


def _offset_disks(radius: float = RADIUS_PX):
    r = int(np.ceil(radius))
    out = []
    for dy in range(-r, r + 1):
        for dx in range(-r, r + 1):
            k = float(kernel(np.hypot(dy, dx), radius))
            if k > 0:
                out.append((dy, dx, k))
    return out


_OFFSETS = _offset_disks()


def masked_dti(pred, truth, mask, valid, *, alpha: float = ALPHA, beta: float = BETA,
               radius: float = RADIUS_PX, mask_can_claim: bool = False) -> dict:
    """Distance-weighted Tversky index with the organiser's catalogue mask.

    ``pred``  prediction surface in [0, 1] (float)
    ``truth`` surrogate ground-truth pixels (bool)
    ``mask``  pixel-exact known-catalogue mask (bool) - excluded from the
              false-positive term, and (``mask_can_claim=False``) also prevented
              from earning true-positive credit
    ``valid`` scored domain (bool)

    Returns ``tp, fp, fn, n_truth, dti, emitted, mean_credit``.
    """
    p = np.asarray(pred, dtype=np.float64)
    g = np.asarray(truth, dtype=bool) & np.asarray(valid, dtype=bool)
    m = np.asarray(mask, dtype=bool)
    p = np.where(np.asarray(valid, dtype=bool), p, 0.0)
    if not mask_can_claim:
        p = np.where(m, 0.0, p)
    n = int(g.sum())
    emitted = int((p > 0).sum())
    if n == 0:
        return dict(tp=0.0, fp=float(p.sum()), fn=0.0, n_truth=0, dti=0.0,
                    emitted=emitted, mean_credit=0.0)
    h, w = p.shape
    yy, xx = np.nonzero(g)
    cred = np.zeros(yy.size, dtype=np.float64)
    for dy, dx, k in _OFFSETS:
        ny, nx = yy + dy, xx + dx
        ok = (ny >= 0) & (ny < h) & (nx >= 0) & (nx < w)
        cred[ok] = np.maximum(cred[ok], p[ny[ok], nx[ok]] * k)
    tp = float(cred.sum())
    dist = distance_transform_edt(~g)
    fp = float((p * (1.0 - kernel(dist, radius))).sum())
    fn = float(n) - tp
    dti = tp / (tp + alpha * fp + beta * fn + EPS)
    return dict(tp=tp, fp=fp, fn=fn, n_truth=n, dti=float(dti), emitted=emitted,
                mean_credit=float(tp / max(emitted, 1)))


def masked_dti_binary(pred_bool, truth, mask, valid, **kw) -> dict:
    return masked_dti(np.asarray(pred_bool, dtype=np.float64), truth, mask, valid, **kw)


def binary_dti_kdtree(dots, truth, mask, valid, *, radius: float = RADIUS_PX) -> dict:
    """Exact masked DTI for a *binary* field, evaluated with KD-trees.

    Mathematically identical to :func:`masked_dti` for a {0,1} field, but O(n)
    in the number of dots/truth pixels instead of a full-grid distance
    transform, so a whole emission curve can be evaluated in seconds.
    """
    from scipy.spatial import cKDTree

    d = np.asarray(dots, bool) & np.asarray(valid, bool) & ~np.asarray(mask, bool)
    g = np.asarray(truth, bool) & np.asarray(valid, bool)
    n = int(g.sum())
    emitted = int(d.sum())
    if n == 0 or emitted == 0:
        return dict(tp=0.0, fp=float(emitted), fn=float(n), n_truth=n, dti=0.0, emitted=emitted,
                    mean_credit=0.0, kappa=0.0)
    dy, dx = np.nonzero(d)
    ty, tx = np.nonzero(g)
    td = cKDTree(np.column_stack([dy, dx]))
    dt = cKDTree(np.column_stack([ty, tx]))
    dtru, _ = td.query(np.column_stack([ty, tx]), k=1)
    tp = float(kernel(dtru, radius).sum())
    ddot, _ = dt.query(np.column_stack([dy, dx]), k=1)
    kappa_vec = kernel(ddot, radius)
    fp = float((1.0 - kappa_vec).sum())
    fn = float(n) - tp
    dti = tp / (tp + ALPHA * fp + BETA * fn + EPS)
    return dict(tp=tp, fp=fp, fn=fn, n_truth=n, dti=float(dti), emitted=emitted,
                mean_credit=float(tp / emitted), kappa=float(kappa_vec.sum() / emitted))


def tp_only_binary(pred_bool, truth, valid) -> float:
    """TP credit of a binary field against a truth set (no FP term)."""
    p = np.asarray(pred_bool, dtype=bool) & np.asarray(valid, dtype=bool)
    g = np.asarray(truth, dtype=bool) & np.asarray(valid, dtype=bool)
    if not g.any():
        return 0.0
    dist = distance_transform_edt(~p)
    return float(kernel(dist[g]).sum())


def sgmc_off_catalogue(sgmc_raster: np.ndarray, catalogue: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Truth = SGMC fault pixels that are not catalogue pixels (inside footprint)."""
    sg = np.asarray(sgmc_raster)
    sg = (sg > 0) & (sg != 255) & (sg != np.iinfo(sg.dtype).max if np.issubdtype(sg.dtype, np.integer) else True)
    return sg & ~np.asarray(catalogue, bool) & np.asarray(valid, bool)


def isolated_catalogue_components(catalogue: np.ndarray, valid: np.ndarray,
                                  min_separation_px: int = 6, min_size_px: int = 12) -> dict:
    """Split the catalogue into (mask, truth) so the mask does not reveal the truth.

    A connected component of catalogue pixels becomes *truth* only if every one
    of its pixels is at least ``min_separation_px`` from every pixel of every
    other component; the remaining catalogue is the mask.  Predicting the mask
    flank is therefore uninformative for the held-out truth.
    """
    from scipy.spatial import cKDTree

    cat = np.asarray(catalogue, bool) & np.asarray(valid, bool)
    lab, n = cc_label(cat, structure=np.ones((3, 3), dtype=int))
    if n == 0:
        raise ValueError("empty catalogue")
    sizes = np.bincount(lab.ravel())
    yy, xx = np.nonzero(cat)
    labels = lab[yy, xx]
    pts = np.column_stack([yy, xx]).astype(np.float64)
    tree = cKDTree(pts)
    # Nearest pixel that does NOT belong to the same component.
    k = min(40, max(8, int(np.ceil(2 * np.sqrt(int(sizes.max()))))))
    dist, idx = tree.query(pts, k=k, workers=-1)
    other_min = np.full(yy.size, np.inf)
    for j in range(k):
        diff = labels[idx[:, j]] != labels
        other_min = np.where(diff & (dist[:, j] < other_min), dist[:, j], other_min)
    isolated_px = other_min >= float(min_separation_px)
    truth = np.zeros_like(cat)
    for i in range(1, n + 1):
        if sizes[i] < min_size_px:
            continue
        sel = labels == i
        if isolated_px[sel].all():
            truth[yy[sel], xx[sel]] = True
    return dict(truth=truth, mask=cat & ~truth,
                n_components=int(n), truth_px=int(truth.sum()),
                mask_px=int((cat & ~truth).sum()))


def random_dots(valid: np.ndarray, mask: np.ndarray, n: int, seed: int = 0,
                min_separation_px: float = 0.0) -> np.ndarray:
    """Uniform random binary control of exactly ``n`` pixels off the mask."""
    rng = np.random.default_rng(seed)
    allowed = np.asarray(valid, bool) & ~np.asarray(mask, bool)
    idx = np.flatnonzero(allowed.ravel())
    n = min(int(n), idx.size)
    pick = rng.choice(idx, size=n, replace=False)
    out = np.zeros(allowed.size, dtype=bool)
    out[pick] = True
    if min_separation_px > 0:
        from .lattice import greedy_min_separation
        s = np.zeros(allowed.size, dtype=np.float64)
        s[pick] = 1.0
        out = greedy_min_separation(s.reshape(allowed.shape), allowed, n_max=n,
                                    min_separation_px=min_separation_px)
    return out.reshape(allowed.shape)


def flank_ring(catalogue: np.ndarray, valid: np.ndarray, inner_px: int = 1,
               outer_px: int = 3, mask_outside: bool = True) -> np.ndarray:
    """Pixels at ``inner_px`` .. ``outer_px`` from the catalogue (a control)."""
    cat = np.asarray(catalogue, bool)
    d = distance_transform_edt(~cat)
    ring = (d > inner_px) & (d <= outer_px) & np.asarray(valid, bool)
    if mask_outside:
        ring &= ~cat
    return ring


def coverage_summary(pred_bool, truth, mask, valid) -> dict:
    """Coverage of the surrogate truth by a binary field, plus redundancy."""
    r = masked_dti_binary(pred_bool, truth, mask, valid)
    emitted = int(r["emitted"])
    if emitted:
        dil = binary_dilation(np.asarray(pred_bool, bool), structure=np.ones((3, 3), bool))
        # a truth pixel is 'covered' if a prediction dot lies within 3 px
        d = distance_transform_edt(~np.asarray(pred_bool, bool))
        cov = float(((np.asarray(truth, bool) & (d <= 3.0)).sum()) / max(int(np.asarray(truth, bool).sum()), 1))
    else:
        dil = pred_bool
        cov = 0.0
    return dict(**r, truth_covered_frac=cov)
