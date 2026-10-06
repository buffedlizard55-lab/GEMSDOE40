"""Spatially blocked holdout instruments used to decide what may be shipped.

Two instruments, both implemented from the metric that is defined in ``gems40.metric``:

1. ``LM`` -- **the decision instrument**.  Spatially blocked (4 quadrants), truth set that is
   off-catalogue *by construction*: the USGS State Geologic Map Compilation (SGMC) fault pixels
   lying more than 300 m from every catalogue pixel, inside an eroded test quadrant.  Because the
   truth excludes the catalogue, mass spent on the catalogue and its flank is charged alpha and
   earns nothing -- which is what the live scorer does to it (organiser clarification:
   https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516 ).
   This instrument and its prevalence calibration were introduced by the sibling repository
   ``GEMSDOE32`` (``src/gems32/live_mirror.py``); the definition is reproduced here so that every
   number in this repository is comparable with that work, and its live record is re-verified
   below before it is used (``validate_instrument``).

2. ``CAT-HID`` -- the *conventional* holdout: 20 % of catalogue fault components are hidden and
   scored in an eroded quadrant with a 1.5 km collar.  It is reported for continuity, but it
   cannot rank live scores (its truth is the catalogue itself), so it is never used to decide.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import rasterio
from scipy.ndimage import binary_dilation, binary_erosion, label

from .metric import ALPHA, BETA, EPS, kernel

FOLD_NAMES = ("NW", "NE", "SW", "SE")
FOLDS = (0, 1, 2, 3)
COLLAR_PX = 15        # 1.5 km buffer between the visible catalogue and the test quadrant
DOMAIN_ERODE = 12     # 1.2 km boundary erosion of the test quadrant
CAT_FLANK_PX = 3      # 300 m: the metric's own kernel radius
G_LIVE_TOTAL = 12691.0  # hidden-truth prevalence px, calibrated by the sibling repo from live scores


def read_binary(path: str | Path) -> np.ndarray:
    with rasterio.open(path) as ds:
        a = ds.read(1)
    return np.isfinite(a) & (a > 0)


def quadrant_ids(footprint: np.ndarray) -> np.ndarray:
    foot = np.asarray(footprint, bool)
    yy, xx = np.nonzero(foot)
    ym, xm = int(np.median(yy)), int(np.median(xx))
    H, W = foot.shape
    gy, gx = np.ogrid[:H, :W]
    q = np.full((H, W), -1, np.int8)
    q[(gy < ym) & (gx < xm) & foot] = 0
    q[(gy < ym) & (gx >= xm) & foot] = 1
    q[(gy >= ym) & (gx < xm) & foot] = 2
    q[(gy >= ym) & (gx >= xm) & foot] = 3
    return q


def _bbox_of(mask: np.ndarray, pad: int = 6) -> tuple[slice, slice]:
    rows = np.flatnonzero(mask.any(axis=1))
    cols = np.flatnonzero(mask.any(axis=0))
    H, W = mask.shape
    return (slice(max(0, int(rows[0]) - pad), min(H, int(rows[-1]) + pad + 1)),
            slice(max(0, int(cols[0]) - pad), min(W, int(cols[-1]) + pad + 1)))


@dataclass
class Cell:
    key: str
    fold: int
    bbox: tuple[slice, slice]
    domain: np.ndarray
    truth: np.ndarray


@dataclass
class LiveMirror:
    foot: np.ndarray
    labels: np.ndarray
    sgmc_off: np.ndarray
    quad: np.ndarray
    cells: list[Cell]


@dataclass
class CatHidden:
    foot: np.ndarray
    labels: np.ndarray
    quad: np.ndarray
    cells: list[Cell]


def load_live_mirror(data_dir: str | Path, sgmc_name: str = "derived_sgmc_faults_100m_u8.tif",
                     external_subdir: str = "external") -> LiveMirror:
    """Build the four blocked off-catalogue cells of the LM instrument."""
    ddir = Path(data_dir)
    sample = ddir / "example_submission.tif"
    if not sample.exists():
        sample = ddir / "sample_submission.tif"
    with rasterio.open(sample) as ds:
        foot = np.isfinite(ds.read(1))
    labels = read_binary(ddir / "existing_faults.tif") & foot
    p = ddir / external_subdir / sgmc_name
    if not p.exists():
        p = ddir / sgmc_name
    with rasterio.open(p) as ds:
        sgmc = ds.read(1) > 0
    sgmc_off = sgmc & foot & ~labels & ~binary_dilation(labels, iterations=CAT_FLANK_PX)
    quad = quadrant_ids(foot)
    cells: list[Cell] = []
    for fold in FOLDS:
        q = quad == fold
        sl = _bbox_of(q)
        domain = binary_erosion(q, iterations=DOMAIN_ERODE)
        cells.append(Cell(key=f"fold{FOLD_NAMES[fold]}", fold=fold, bbox=sl,
                          domain=domain[sl], truth=(sgmc_off & domain)[sl]))
    return LiveMirror(foot=foot, labels=labels, sgmc_off=sgmc_off, quad=quad, cells=cells)


def load_cat_hidden(data_dir: str | Path, hide_frac: float = 0.20, seed: int = 21) -> CatHidden:
    """Build the four blocked catalogue-component holdout cells (conventional proxy)."""
    ddir = Path(data_dir)
    sample = ddir / "example_submission.tif"
    if not sample.exists():
        sample = ddir / "sample_submission.tif"
    with rasterio.open(sample) as ds:
        foot = np.isfinite(ds.read(1))
    labels = read_binary(ddir / "existing_faults.tif") & foot
    quad = quadrant_ids(foot)
    comp, n_comp = label(labels, structure=np.ones((3, 3), int))
    rng = np.random.default_rng(seed)
    ids = np.arange(1, n_comp + 1)
    sizes = np.bincount(comp.ravel(), minlength=n_comp + 1)
    perm = rng.permutation(ids)
    cum = np.cumsum(sizes[perm])
    k = int(np.searchsorted(cum, hide_frac * sizes.sum())) + 1
    hidden_ids = perm[:min(k, perm.size)]
    hidden = np.isin(comp, hidden_ids) & foot
    cells: list[Cell] = []
    for fold in FOLDS:
        q = quad == fold
        sl = _bbox_of(q)
        collar = binary_dilation(q, structure=np.ones((3, 3), bool), iterations=COLLAR_PX) & foot
        visible = labels & (q | collar)
        domain = binary_erosion(q, iterations=DOMAIN_ERODE)
        truth = hidden & q & domain & ~binary_dilation(visible, iterations=CAT_FLANK_PX)
        cells.append(Cell(key=f"fold{FOLD_NAMES[fold]}", fold=fold, bbox=sl,
                          domain=domain[sl], truth=truth[sl]))
    return CatHidden(foot=foot, labels=labels, quad=quad, cells=cells)


def lm_score(mask: np.ndarray, ctx: LiveMirror, calibrate_prevalence: bool = True,
             g_lb_total: float = G_LIVE_TOTAL) -> dict:
    """Score a binary support mask under LM (prevalence-calibrated by default)."""
    m = np.asarray(mask, bool) & ctx.foot
    foot_px = float(ctx.foot.sum())
    per_fold, per_fold_cal, detail = {}, {}, {}
    for cell in ctx.cells:
        sl = cell.bbox
        p = m[sl] & cell.domain
        g = cell.truth
        n_p, n_g = int(p.sum()), int(g.sum())
        if n_p == 0 or n_g == 0:
            per_fold[cell.key] = per_fold_cal[cell.key] = 0.0
            detail[cell.key] = dict(n_p=n_p, n_g=n_g)
            continue
        from scipy.ndimage import distance_transform_edt
        dp = distance_transform_edt(~p)
        dg = distance_transform_edt(~g)
        tp_p = float(kernel(dp[g]).sum())     # credit each truth pixel gets from the nearest dot
        tp_g = float(kernel(dg[p]).sum())     # kernel mass each dot delivers
        tp_w = 0.5 * (tp_p + tp_g)
        denom = ALPHA * n_p + BETA * n_g + BETA * (tp_g - tp_p) + EPS
        per_fold[cell.key] = tp_w / denom
        if calibrate_prevalence:
            g_cal = g_lb_total * (float(cell.domain.sum()) / foot_px)
            denom_cal = ALPHA * n_p + BETA * g_cal + BETA * (tp_g - tp_p) + EPS
            per_fold_cal[cell.key] = tp_w / denom_cal
        else:
            per_fold_cal[cell.key] = per_fold[cell.key]
        detail[cell.key] = dict(n_p=n_p, n_g=n_g, tp_w=tp_w)
    return dict(lm_mean=float(np.mean(list(per_fold.values()))),
                lm_calibrated_mean=float(np.mean(list(per_fold_cal.values()))),
                lm_per_fold=per_fold, lm_calibrated_per_fold=per_fold_cal,
                emitted_pixels=int(m.sum()),
                on_catalogue_pixels=int((m & ctx.labels).sum()),
                fold_detail=detail)


def cat_hidden_score(mask: np.ndarray, ctx: CatHidden) -> dict:
    """Score a binary support mask under the conventional catalogue-component holdout."""
    m = np.asarray(mask, bool) & ctx.foot
    from scipy.ndimage import distance_transform_edt
    vals = {}
    for cell in ctx.cells:
        sl = cell.bbox
        p = m[sl] & cell.domain
        g = cell.truth
        n_p, n_g = int(p.sum()), int(g.sum())
        if n_p == 0 or n_g == 0:
            vals[cell.key] = 0.0
            continue
        dp = distance_transform_edt(~p)
        dg = distance_transform_edt(~g)
        tp_p = float(kernel(dp[g]).sum())
        tp_g = float(kernel(dg[p]).sum())
        tp_w = 0.5 * (tp_p + tp_g)
        vals[cell.key] = tp_w / (ALPHA * n_p + BETA * n_g + BETA * (tp_g - tp_p) + EPS)
    return dict(cat_hidden_mean=float(np.mean(list(vals.values()))), cat_hidden_per_fold=vals,
                emitted_pixels=int(m.sum()))
