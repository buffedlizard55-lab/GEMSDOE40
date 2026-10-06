"""Metric-optimal emission for the GEMS distance-weighted Tversky index.

Why this module exists
----------------------
The organizer's metric (problem description, page 967) is

    DTI = TPw / (TPw + alpha*FPw + beta*FNw + eps),  alpha = 0.2, beta = 0.8

with the triangular kernel k(d) = max(1 - d/R, 0), R = 300 m = 3 px.

Two algebraic facts follow, and both are used here (and verified against the
independent implementation in :mod:`gemsdoe40.research_metric` by the tests):

1. **Binarization weakly dominates.**  Fix the support of a prediction and scale
   every value by lambda <= 1.  Then TPw -> lambda*TPw and FPw -> lambda*FPw, so

       DTI(lambda) = lambda*T / (0.8*N + 0.2*lambda*T + 0.2*lambda*F)

   whose derivative in lambda is 0.8*N*T / (denominator^2) > 0.  For a fixed
   support the score is therefore maximized at lambda = 1, i.e. by writing
   exactly 1.0 on every emitted pixel.  A graded confidence raster is
   *weakly dominated* by its own binarization.

2. **For a binary emission the score has an exact closed form.**  With p in
   {0,1}, M = sum(p) = number of emitted pixels, T = TPw, and

       C_pred = sum_x p(x) * max_g k(d(x,g))     ("credit seen from the dots")

   the false-positive term is FPw = M - C_pred, because every emitted pixel
   contributes its own triangular kernel value to the penalty *unless* it sits
   on a ground-truth pixel.  Substituting into the definition gives

       DTI = T / (0.8*N + 0.2*M + 0.2*(T - C_pred))                     (★)

   with N = number of hidden test pixels.  (★) is verified to 1e-9 against the
   independent implementation of the published equations in
   ``gemsdoe40.research_metric`` by ``tests/test_emission.py``.  Two special
   cases: a dot set with one dot per covered truth pixel has C_pred = T and the
   correction vanishes; a uniform dense emission has C_pred >> T and the
   correction lowers the score, which is the algebraic reason a blind lattice is
   so inefficient.  The score is *credit per unit mass* against a fixed budget;
   arrangement matters only through T and C_pred.

   ``binary_dti_thin`` returns (★) evaluated at C_pred = T.  It is an
   approximation and is labelled as one everywhere it is used.  It uses the
   pixel-exact exclusion of the known-fault mask, which the organizer confirmed
   applies to the penalty terms and which the sibling project measured directly:
   adding the whole catalogue to an off-catalogue emission left the live score
   unchanged at 0.1563.

3. **The credit bar.**  Differentiating (★) for one extra emitted pixel whose
   expected credit to the nearest hidden pixel is k, and whose mass cost is 1:

       d DT I > 0  <=>  k > 0.2 * DTI

   This is the exact break-even rule used to stop emitting.  (GEMSDOE32 derived
   the same bar independently as "credit bar k > 0.2*DTI"; at its live 0.26 that
   is 0.052, matching its empirically measured 0.0548.)

Consequences used by the generator
----------------------------------
* Emit a *binary* field.
* Rank candidate pixels by expected credit and stop at the bar.
* Thin aggressively: two dots within R of each other cannot both raise T, but
  both cost mass.

Nothing in this module is a score.  ``model_dti`` is a model.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

ALPHA = 0.2
BETA = 0.8
RADIUS_PX = 3.0


def binary_dti_exact(
    total_credit: float,
    pred_credit: float,
    mass: int | float,
    n_truth: float,
) -> float:
    """Exact metric value for a binary, off-mask emission.  Implements (★).

    ``total_credit`` is TPw (credit accumulated at ground-truth pixels),
    ``pred_credit`` is ``sum_x p(x) * max_g k(d(x,g))`` (credit seen from the
    emitted pixels), ``mass`` is the number of emitted pixels, ``n_truth`` is the
    number of ground-truth pixels in the scored subset.
    """
    if mass < 0 or n_truth < 0 or total_credit < 0 or pred_credit < 0:
        raise ValueError("all arguments must be non-negative")
    denom = BETA * float(n_truth) + ALPHA * float(mass) + ALPHA * (float(total_credit) - float(pred_credit))
    if denom <= 0:
        return 0.0
    return float(total_credit) / denom


def binary_dti_thin(total_credit: float, mass: int | float, n_truth: float) -> float:
    """Exact (★) evaluated at ``C_pred = T`` -- the one-dot-per-truth-pixel case.

    ``TPw / (0.8*N + 0.2*M)``.  Documented as an approximation; the generator
    measures the correction term on the proxy truth before using it.
    """
    return binary_dti_exact(total_credit, total_credit, mass, n_truth)


def pred_credit_of_binary(prediction: np.ndarray, truth: np.ndarray,
                          *, radius_px: float = RADIUS_PX) -> float:
    """``C_pred = sum_x p(x) * max_g k(d(x,g))`` via a Euclidean distance transform."""
    p = np.asarray(prediction, dtype=np.float32)
    g = np.asarray(truth, dtype=bool)
    if p.shape != g.shape:
        raise ValueError("prediction and truth must share a shape")
    if not g.any():
        return 0.0
    from scipy import ndimage

    d = ndimage.distance_transform_edt(~g)
    k = np.maximum(1.0 - d / float(radius_px), 0.0).astype(np.float32)
    return float((p * k).sum(dtype=np.float64))


def credit_bar(dti: float) -> float:
    """Break-even expected credit for one more emitted pixel: ``0.2 * DTI``."""
    return ALPHA * float(dti)


def metric_required_credit(dti: float, mass: int | float, n_truth: float) -> float:
    """Total credit ``T`` needed to reach ``dti`` at this mass and truth count.

    From (★):  T = dti * (0.8*N + 0.2*M).
    """
    return float(dti) * (BETA * float(n_truth) + ALPHA * float(mass))


def metric_max_mass(credit: float, n_truth: float, dti: float) -> float:
    """Largest mass that reaches ``dti`` with a fixed total credit.

    From (★) at C_pred = T:  M = (T / dti - 0.8*N) / 0.2.  Negative when the target
    is out of reach at that credit.
    """
    if dti <= 0:
        raise ValueError("dti must be positive")
    return (float(credit) / float(dti) - BETA * float(n_truth)) / ALPHA


@dataclass
class Emission:
    """A binary emission and the evidence measured for it."""

    values: np.ndarray          # float32, exactly 0.0 / 1.0 (NaN never inside support)
    spacing_px: float
    mass: int
    rejected_by_flank: int = 0

    def __post_init__(self) -> None:
        if self.values.dtype != np.float32:
            raise TypeError("emission values must be float32")


def value_ranked_thinning(
    score: np.ndarray,
    support: np.ndarray,
    *,
    spacing_px: float,
    max_mass: int | None = None,
    min_score: float = 0.0,
) -> np.ndarray:
    """Greedy value-ranked max-separation (Poisson-disk) selection.

    Candidates inside ``support`` with ``score > min_score`` are visited in
    descending score order.  A candidate is accepted when no already-accepted
    pixel lies within ``spacing_px``; accepting a pixel excludes every pixel of
    its disk of that radius.  This is the standard greedy solution of the
    weighted max-separation packing problem, and it is what makes the emission
    mass-efficient: two dots inside the metric's own radius R cannot both raise
    TP, but both are charged to FP.

    Returns a float32 array with 1.0 on accepted pixels and 0.0 elsewhere.
    """
    if spacing_px <= 0:
        raise ValueError("spacing_px must be positive")
    s = np.asarray(score, dtype=np.float64)
    sup = np.asarray(support, dtype=bool)
    if s.shape != sup.shape:
        raise ValueError("score and support must share a shape")
    h, w = s.shape
    cand = sup & np.isfinite(s) & (s > min_score)
    n_cand = int(cand.sum())
    out = np.zeros((h, w), dtype=np.float32)
    if n_cand == 0:
        return out

    ridx, cidx = np.nonzero(cand)
    vals = s[ridx, cidx]
    order = np.argsort(-vals, kind="stable")
    ridx = ridx[order]
    cidx = cidx[order]

    r = float(spacing_px)
    rad = int(np.ceil(r))
    blocked = np.zeros((h, w), dtype=bool)
    offs = [(dy, dx) for dy in range(-rad, rad + 1) for dx in range(-rad, rad + 1)
            if dy * dy + dx * dx <= r * r]
    accepted = 0
    for i in range(ridx.size):
        y = int(ridx[i])
        x = int(cidx[i])
        if blocked[y, x]:
            continue
        out[y, x] = 1.0
        accepted += 1
        for dy, dx in offs:
            yy = y + dy
            xx = x + dx
            if 0 <= yy < h and 0 <= xx < w:
                blocked[yy, xx] = True
        if max_mass is not None and accepted >= max_mass:
            break
    return out


# --------------------------------------------------------------------------- #
# live-score-calibrated instrument (reproduced from the sibling repository)
# --------------------------------------------------------------------------- #

#: Saturation constant and surrogate-to-real credit mapping of the published
#: instrument, least-squares fitted by ``GEMSDOE39`` (session 2,
#: ``registry/h40_report_h40e-30k.json``, ``src/gems39/h40.py``) to the 12
#: organizer-scored artifacts whose rasters exist in this programme's prior
#: corpus: RMSE 0.0364, in-sample Spearman +0.923, leave-one-out Spearman
#: +0.902, leave-one-out RMSE 0.0397.  The fit is *not* ours; it is reproduced
#: here so every number in this repository can be audited against the same
#: live-anchored model, and the ``w`` it consumes is defined identically.
INSTRUMENT_K = 5796.844642382955
INSTRUMENT_A = 1.1689948411387971


def dot_credit_w(dots: np.ndarray, truth: np.ndarray, radius_px: float = 3.0) -> float:
    """Mean kernel credit per emitted dot against ``truth`` (the instrument's ``w``).

    ``w = mean over emitted dots of max(0, 1 - d/R)`` with ``d`` the distance to
    the nearest truth pixel and ``R = 300 m = 3 px``.  This is a *dot-centric*
    quantity (the mirror of :func:`credit_components`' truth-centric ``tp``):
    ``w * mass`` is the prediction-side credit ``C_pred``.  Verified in this
    repository against the nine published ``w`` values of the sibling fit
    (max absolute difference < 1e-4 on the files present locally).
    """
    from scipy.ndimage import distance_transform_edt

    g = np.asarray(truth, dtype=bool)
    d = np.asarray(dots) > 0
    n = int(d.sum())
    if n == 0 or not g.any():
        return 0.0
    k = np.maximum(1.0 - distance_transform_edt(~g) / float(radius_px), 0.0)
    return float(k[d].mean())


def calibrated_live_score(n_dots: int | float, w: float,
                          K: float = INSTRUMENT_K, a: float = INSTRUMENT_A) -> float:
    """Predicted organizer score from (dot count, per-dot surrogate credit).

    Saturating-coverage model fitted to 12 organizer-scored artifacts::

        TP_hat = K * (1 - exp(-a * n * w / K))
        S_hat  = TP_hat / (0.2 * TP_hat + 0.2 * n + 0.8 * K)

    The saturation is what reproduces the measured negative rank correlation
    between emission size and live score (rho = -0.867 across the 12 anchors),
    which a constant-credit model cannot express.  ``w`` is measured on the
    off-catalogue surrogate truth (SGMC faults outside the provided catalogue);
    the model is fitted on ``w`` in [0.048, 0.105], so values far outside that
    range are extrapolations and are flagged wherever they are reported.

    Break-even corollary (their budget rule, reproduced): a marginal dot pays iff
    ``a * exp(-a * C_sur / K) * w_marginal > 0.2 * S_hat``.
    """
    n = float(n_dots)
    tp = K * (1.0 - np.exp(-a * n * float(w) / K))
    return float(tp / (0.2 * tp + 0.2 * n + 0.8 * K))


def instrument_break_even(n_dots: int | float, w: float,
                          K: float = INSTRUMENT_K, a: float = INSTRUMENT_A) -> float:
    """Marginal ``w`` a new dot must beat for the instrument score to rise.

    Exact marginal condition of the closed form: with
    ``dTP/dn = a * exp(-a*n*w/K) * w_marginal`` and denominator
    ``den = 0.2*TP + 0.2*n + 0.8*K``,
    ``dS/dn > 0  <=>  dTP * (0.2*n + 0.8*K) > 0.2 * TP``.
    The sibling repository states the approximation
    ``a*exp(-a*n*w/K)*w_marginal > 0.2*S``, which differs by the factor
    ``den / (0.2*n + 0.8*K)`` (1.0 only when ``TP = 0``); the test suite pins the
    exact form.
    """
    n = float(n_dots)
    tp = K * (1.0 - np.exp(-a * n * float(w) / K))
    d_tp_per_w = a * np.exp(-a * n * float(w) / K)
    return float(0.2 * tp / ((0.2 * n + 0.8 * K) * d_tp_per_w))


def thinning_curve_snapshots(
    score: np.ndarray,
    support: np.ndarray,
    *,
    spacing_px: float,
    breakpoints: list[int],
    on_snapshot,
    min_score: float = 0.0,
) -> None:
    """One streaming pass of the value-ranked thinning, snapshotted at masses.

    The greedy selection is nested in the mass budget: the set accepted for
    budget ``m1`` is a subset of the set accepted for ``m2 > m1``.  Re-running
    the selection once per budget therefore wastes work, and this routine makes
    the u(M) curve affordable in one pass.  ``on_snapshot(mass, values)`` is
    called with a fresh float32 array at every breakpoint.
    """
    if spacing_px <= 0:
        raise ValueError("spacing_px must be positive")
    bp = sorted(int(b) for b in breakpoints)
    if not bp or bp[0] <= 0:
        raise ValueError("breakpoints must be positive")
    s = np.asarray(score, dtype=np.float64)
    sup = np.asarray(support, dtype=bool)
    if s.shape != sup.shape:
        raise ValueError("score and support must share a shape")
    h, w = s.shape
    cand = sup & np.isfinite(s) & (s > min_score)
    if not cand.any():
        return
    ridx, cidx = np.nonzero(cand)
    order = np.argsort(-s[ridx, cidx], kind="stable")
    ridx = ridx[order]
    cidx = cidx[order]
    r = float(spacing_px)
    rad = int(np.ceil(r))
    offs = [(dy, dx) for dy in range(-rad, rad + 1) for dx in range(-rad, rad + 1)
            if dy * dy + dx * dx <= r * r]
    blocked = np.zeros((h, w), dtype=bool)
    out = np.zeros((h, w), dtype=np.float32)
    accepted = 0
    nxt = 0
    for i in range(ridx.size):
        y = int(ridx[i])
        x = int(cidx[i])
        if blocked[y, x]:
            continue
        out[y, x] = 1.0
        accepted += 1
        for dy, dx in offs:
            yy = y + dy
            xx = x + dx
            if 0 <= yy < h and 0 <= xx < w:
                blocked[yy, xx] = True
        while nxt < len(bp) and accepted >= bp[nxt]:
            on_snapshot(accepted, out.copy())
            nxt += 1
        if nxt >= len(bp):
            return
    while nxt < len(bp):
        on_snapshot(accepted, out.copy())
        nxt += 1


def credit_of_binary(
    prediction: np.ndarray,
    truth: np.ndarray,
    *,
    radius_px: float = RADIUS_PX,
    mask: np.ndarray | None = None,
) -> float:
    """TPw of a binary prediction under the organizer's triangular kernel.

    Implemented by brute force over the (2R+1)^2 neighbourhood, using an
    offset-shift maximum -- an independent route from the distance-transform
    implementation in :mod:`gemsdoe40.research_metric`, so agreement between the
    two is a real check.
    """
    p = np.asarray(prediction, dtype=np.float32)
    g = np.asarray(truth, dtype=bool)
    if p.shape != g.shape:
        raise ValueError("prediction and truth must share a shape")
    if mask is not None:
        m = np.asarray(mask, dtype=bool)
        if m.shape != p.shape:
            raise ValueError("mask must match the prediction shape")
        p = np.where(m, 0.0, p)
        g = g & ~m
    if not g.any():
        return 0.0
    h, w = p.shape
    rad = int(np.floor(radius_px))
    best = np.zeros((h, w), dtype=np.float32)
    for dy in range(-rad, rad + 1):
        for dx in range(-rad, rad + 1):
            d = float(np.hypot(dx, dy))
            if d > radius_px:
                continue
            wt = np.float32(1.0 - d / radius_px)
            shifted = np.zeros_like(p)
            ys = slice(max(0, dy), min(h, h + dy))
            yd = slice(max(0, -dy), min(h, h - dy))
            xs = slice(max(0, dx), min(w, w + dx))
            xd = slice(max(0, -dx), min(w, w - dx))
            shifted[ys, xs] = p[yd, xd]
            np.maximum(best, shifted * wt, out=best)
    return float(best[g].sum(dtype=np.float64))


def credit_density_curve(
    score: np.ndarray,
    support: np.ndarray,
    truth: np.ndarray,
    *,
    spacing_px: float,
    mass_grid: list[int],
    mask: np.ndarray | None = None,
) -> list[dict]:
    """Credit density u(M) = T/M along the value-ranked thinning curve.

    ``mass_grid`` must be ascending.  Each entry thins the same ranking up to
    the requested mass, so the curve is monotone in M and directly comparable
    across emissions built from the same field.
    """
    m_prev = 0
    out: list[dict] = []
    for m in mass_grid:
        if m <= m_prev:
            raise ValueError("mass_grid must be strictly ascending")
        emit = value_ranked_thinning(score, support, spacing_px=spacing_px, max_mass=m)
        mass = int(emit.sum())
        credit = credit_of_binary(emit, truth, mask=mask)
        pred_credit = pred_credit_of_binary(emit, truth)
        n_truth = int(truth.sum()) if mask is None else int((truth & ~mask).sum())
        out.append({
            "spacing_px": float(spacing_px),
            "mass": mass,
            "credit": credit,
            "pred_credit": pred_credit,
            "n_truth": n_truth,
            "u": credit / mass if mass else 0.0,
            "correction": ALPHA * (credit - pred_credit),
            "dti_proxy_exact": binary_dti_exact(credit, pred_credit, mass, n_truth),
        })
        m_prev = m
    return out


def model_dti(
    curves: list[dict],
    *,
    n_truth: float,
    calibration: float,
    truth_label: str,
) -> list[dict]:
    """Apply the live-anchored calibration to a credit-density curve.

    ``calibration`` is u_live / u_proxy measured on artifacts that already have
    a reported live score; it converts proxy credit density into expected live
    credit density.  The returned ``dti`` is a MODEL, never a score.
    """
    out: list[dict] = []
    for row in curves:
        u_live = float(row["u"]) * float(calibration)
        credit = u_live * float(row["mass"])
        pred_credit = float(row.get("pred_credit", row["credit"])) * float(calibration)
        out.append({
            **row,
            "truth": truth_label,
            "calibration": float(calibration),
            "n_truth": float(n_truth),
            "u_live_model": u_live,
            "credit_live_model": credit,
            "dti_model": binary_dti_exact(credit, pred_credit, row["mass"], n_truth),
            "dti_model_thin": binary_dti_thin(credit, row["mass"], n_truth),
        })
    return out


def maximin_choice(
    rows: list[dict],
    *,
    key: str = "dti_model",
    group_keys: tuple[str, ...] = (),
) -> dict:
    """Pick the emission that maximizes the worst case over model parameters.

    ``group_keys`` names the columns that identify one emission geometry (e.g.
    spacing); every scenario for that geometry (different N and calibration) is
    then reduced by its minimum, and the geometry with the largest minimum is
    chosen.  Maximin is the right rule here because the parameters (hidden-truth
    count, proxy calibration) are uncertain and a submission slot is scarce.
    """
    groups: dict[tuple, list[float]] = {}
    meta: dict[tuple, dict] = {}
    for r in rows:
        g = tuple(r[k] for k in group_keys) if group_keys else ("all",)
        groups.setdefault(g, []).append(float(r[key]))
        meta.setdefault(g, r)
    best = None
    for g, vals in groups.items():
        worst = min(vals)
        if best is None or worst > best["worst_case"]:
            best = {"group": g, "worst_case": worst, "n_scenarios": len(vals),
                    "mean": float(np.mean(vals)), "row": meta[g]}
    if best is None:
        raise ValueError("no rows to choose from")
    return best
