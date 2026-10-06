#!/usr/bin/env python3
"""Euler deconvolution for the H45 depth-clustering candidate.

Implements Reid, Allsop, Granser, Millett & Somerton (1990), *Magnetic interpretation in
three dimensions using Euler deconvolution*, GEOPHYSICS 55(1):80-91,
https://doi.org/10.1190/1.1442774 (author PDF:
https://www.reid-geophys.co.uk/wp-content/uploads/2017/11/Reid-et-al-1990.pdf)
with the structural-index corrections of Reid & Thurston (2014), *The structural index in
gravity and magnetic interpretation: Errors, uses, and abuses*, GEOPHYSICS 79(4):J61-J77,
https://doi.org/10.1190/geo2013-0235.1 (author PDF:
https://www.reid-geophys.co.uk/wp-content/uploads/2017/11/Reid-Thurston-2014.pdf).

Equations, verbatim from the sources
------------------------------------
Reid et al. (1990) eq. (1), structural index N != 0::

    (x-x0) dT/dx + (y-y0) dT/dy + (z-z0) dT/dz = N (B - T)

Reid et al. (1990) eq. (2), structural index N = 0 (contact); the unknown constant ``A``
absorbs amplitude, strike and dip factors::

    (x-x0) dT/dx + (y-y0) dT/dy + (z-z0) dT/dz = A

Reid & Thurston (2014) Table 1, indices used here:

    ======================  =========  =========
    source                  SI (mag)   SI (grav)
    ======================  =========  =========
    thin sheet edge             1          0
    finite contact / fault      0         -1
    ======================  =========  =========

Reid & Thurston state that the gravity SI for a *finite* step is -1 and that it "requires
a more generalized formulation"; the -1 family is therefore run with a disclosed
omitted-variable caveat, and the neighbouring values 0 and 0.5 are run as robustness
families.  Reid et al. (1990) also advise solving a *range* of indices and selecting per
feature by solution clustering, which is what the multi-SI design below does.

Conventions
-----------
* x = east, y = north, z = **down** positive.  Grid rows decrease northward
  (north-up raster with a negative ``dy``), which is handled explicitly.
* The observation plane is z = 0, so a source below the plane has ``z0 > 0``.
* Window coordinates are centred on the window; the solved ``(x0, y0)`` is therefore an
  offset **from the window centre**, in metres.
* ``d/dz`` is the downward-continuation derivative (spectral ``|k|``).

Differences from the H4 solver already in this repository
---------------------------------------------------------
1. Structural indices follow Reid & Thurston (2014) Table 1 instead of a single SI=0 for
   every family: gravity runs at the corrected finite-contact SI of -1 (plus 0, 0.5) and
   magnetics at 0, 0.5 and 1.0, as Reid et al. (1990) explicitly recommend.
2. The unidentifiable along-strike coordinate is handled by an **anisotropic ridge** that
   penalises displacement *along* the local strike only, instead of a hard eigenvalue cut
   with a jump to window-centre anchoring.  Across-strike position stays unpenalised, and
   the penalty fades to zero where the horizontal gradient is isotropic (no strike exists).
3. Small windows (500-900 m) are included so that *near-surface* contacts are resolvable,
   and the cloud weighting is explicitly shallow-selective.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.ndimage import distance_transform_edt, uniform_filter

__all__ = ["Family", "SolveConfig", "SolutionCloud",
           "spectral_derivatives", "solve_euler_grid", "run_family"]


@dataclass(frozen=True)
class Family:
    """One (field, structural index, window) Euler family."""

    name: str                 # band name in training_features.tif
    kind: str                 # "magnetic" | "gravity"
    si: float
    window: int               # window side in grid cells (odd)
    stride: int = 3
    upward_m: float = 0.0     # upward continuation applied before differentiation
    derivative_order: int = 0 # extra vertical derivatives applied to the field


@dataclass(frozen=True)
class SolveConfig:
    dx: float = 100.0
    # A straight contact does not identify the along-strike coordinate: the field barely
    # varies in that direction.  Rather than a hard rank cut or a metre-denominated ridge
    # (whose strength would depend on the field's units), the along-strike constraining
    # power is raised to at least ``min_along_strike_power`` times the across-strike power.
    # With A_along = kappa * A_across the along-strike standard error becomes
    # sigma_across / sqrt(kappa), so kappa = 0.25 means "twice the across-strike slack".
    min_along_strike_power: float = 0.25
    min_depth_m: float = 0.0
    max_depth_m: float = 6000.0
    max_rel_residual: float = 0.30
    max_depth_se_frac: float = 0.35      # conditional depth std err / depth
    max_abs_offset_frac: float = 0.60    # |horizontal offset| / (window/2)
    min_window_valid_frac: float = 0.90
    edge_suppress_cells: int = 6
    # eigenvalue-ratio ramp: ratio <= lo -> full regularisation, ratio >= hi -> none
    isotropic_fade: tuple[float, float] = (0.10, 0.25)


@dataclass
class SolutionCloud:
    """Depth-labelled Euler solutions for one family (all arrays are 1-D and aligned)."""

    family: str
    kind: str
    si: float
    window: int
    col: np.ndarray            # float column (sub-pixel) in the output grid
    row: np.ndarray            # float row (sub-pixel)
    depth: np.ndarray          # metres below the observation plane
    depth_se: np.ndarray       # conditional standard error, metres
    rel_residual: np.ndarray
    offset: np.ndarray         # |horizontal offset from window centre|, metres
    along_strike: np.ndarray   # signed displacement along the local strike, metres
    cross_strike: np.ndarray   # signed displacement across the local strike, metres
    strike_ratio: np.ndarray   # minor/major eigenvalue ratio of the gradient moment
    n_obs: np.ndarray
    meta: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return int(self.col.size)


# --------------------------------------------------------------------------------------
# spectral derivatives
# --------------------------------------------------------------------------------------

def _fill_invalid(values: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Nearest-value fill of invalid cells so the FFT is not polluted by the sentinel.

    The fill is used only to evaluate derivatives; every window overlapping an originally
    invalid cell is rejected later by ``min_window_valid_frac``.
    """
    out = np.where(valid, values, 0.0).astype(np.float64)
    if valid.all():
        return out
    idx = distance_transform_edt(~valid, return_distances=False, return_indices=True)
    return out[tuple(idx)]


def spectral_derivatives(values: np.ndarray, valid: np.ndarray, dx: float,
                         upward_m: float = 0.0, order: int = 0,
                         pad: int = 96) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(F, dF/dx east, dF/dy north, dF/dz down)`` computed consistently.

    All four outputs pass through the *same* spectral operator
    ``exp(-|k| h) * |k|^order``, so the Euler equation is applied to one well-defined
    field and its derivatives of that same field.  A reflect pad suppresses wrap-around.
    """
    f = _fill_invalid(values, valid)
    f = f - f[valid].mean()
    ny, nx = f.shape
    fp = np.pad(f, pad, mode="reflect")
    py, px = fp.shape
    ky = 2.0 * np.pi * np.fft.fftfreq(py, d=dx)[:, None]
    kx = 2.0 * np.pi * np.fft.fftfreq(px, d=dx)[None, :]
    k = np.hypot(ky, kx)
    spec = np.fft.fft2(fp)
    if upward_m:
        spec = spec * np.exp(-k * upward_m)
    base = spec * (k ** order)
    fld = np.real(np.fft.ifft2(base))[pad:pad + ny, pad:pad + nx]
    dfdx = np.real(np.fft.ifft2(base * 1j * kx))[pad:pad + ny, pad:pad + nx]
    # rows decrease northward: d/dy(north) = -d/d(row index)
    dfdy = -np.real(np.fft.ifft2(base * 1j * ky))[pad:pad + ny, pad:pad + nx]
    dfdz = np.real(np.fft.ifft2(base * k))[pad:pad + ny, pad:pad + nx]
    return fld, dfdx, dfdy, dfdz


# --------------------------------------------------------------------------------------
# windowed least squares, fully vectorised with box filters
# --------------------------------------------------------------------------------------

def _box(a: np.ndarray, size: int) -> np.ndarray:
    """Windowed sum over a ``size x size`` centred box, zero padding outside the grid."""
    return uniform_filter(a, size=size, mode="constant", cval=0.0) * float(size * size)


def solve_euler_grid(fx: np.ndarray, fy: np.ndarray, fz: np.ndarray, fval: np.ndarray,
                     si: float, window: int, stride: int, cfg: SolveConfig,
                     n_valid: np.ndarray, rows_per_slab: int = 192) -> dict:
    """Solve Euler's equation for every window on a stride lattice.

    Window sums are accumulated with box filters, so the cost is O(N) in the number of
    grid cells and independent of the window size.

    The grid is processed in horizontal slabs with a halo of ``window`` rows so that peak
    memory stays bounded (the 4x4 normal matrix is 16 floats per cell).  Slab edges use
    zero padding, but the halo rows are discarded, so the result is identical to a
    single-pass computation over the whole grid.

    Coordinate definitions inside a window whose centre is at ``(rw, cw)``:
        ``u = dx * (col - cw)``   east local coordinate (metres)
        ``v = -dx * (row - rw)``  north local coordinate (metres)
    and the Euler right-hand side is ``r = u*fx + v*fy + si*fval``
    (the ``si*fval`` term is absent in the N=0 offset form, where the constant column
    absorbs it).
    """
    ny, nx = fx.shape
    half = window // 2
    row_tops = list(range(half, ny - half, stride))
    if not row_tops:
        raise ValueError(f"window {window} is too large for a {ny}x{nx} grid")
    slabs = []
    start = 0
    while start < len(row_tops):
        chunk = row_tops[start:start + rows_per_slab]
        y0 = max(0, chunk[0] - half)
        y1 = min(ny, chunk[-1] + half + 1)
        slabs.append((chunk, y0, y1))
        start += rows_per_slab

    parts = []
    for chunk, y0, y1 in slabs:
        sl = slice(y0, y1)
        keep = np.searchsorted(np.arange(y0 + half, y1 - half), chunk)
        parts.append(_solve_euler_slab(fx[sl], fy[sl], fz[sl], fval[sl], si, window, stride,
                                       cfg, n_valid[sl], chunk, keep, y0))
    keys = parts[0].keys()
    out = {}
    for k in keys:
        if k in ("row_idx", "col_idx"):
            out[k] = parts[0][k] if k == "col_idx" else np.concatenate([pp[k] for pp in parts])
        else:
            out[k] = np.concatenate([pp[k] for pp in parts], axis=0)
    return out


def _solve_euler_slab(fx: np.ndarray, fy: np.ndarray, fz: np.ndarray, fval: np.ndarray,
                      si: float, window: int, stride: int, cfg: SolveConfig,
                      n_valid: np.ndarray, rows_tops: list, keep: np.ndarray,
                      row_offset: int) -> dict:
    """Single-slab worker for :func:`solve_euler_grid` (see its docstring)."""
    ny, nx = fx.shape
    dx = cfg.dx
    offset_form = abs(si) < 1e-9
    last = 1.0 if offset_form else si

    c_grid = np.broadcast_to(np.arange(nx, dtype=np.float64)[None, :], (ny, nx))
    r_grid = np.broadcast_to(np.arange(ny, dtype=np.float64)[:, None], (ny, nx))
    c2 = c_grid * c_grid
    r2 = r_grid * r_grid
    cr = c_grid * r_grid

    def su(arr):   # sum over window of (col - cw) * arr
        return _box(c_grid * arr, window) - c_grid * _box(arr, window)

    def sr(arr):   # sum over window of (row - rw) * arr
        return _box(r_grid * arr, window) - r_grid * _box(arr, window)

    # design columns
    g1, g2, g3 = fx, fy, fz
    g4 = np.full_like(fx, last)
    design = (g1, g2, g3, g4)

    # ---- normal matrix A[i,j] = sum gi*gj ---------------------------------------------
    a = np.empty((ny, nx, 4, 4), dtype=np.float64)
    for i in range(4):
        for j in range(i, 4):
            a[..., i, j] = _box(design[i] * design[j], window)
            a[..., j, i] = a[..., i, j]

    # ---- right-hand side ---------------------------------------------------------------
    # r = u*fx + v*fy (+ si*fval)
    sb = np.empty((ny, nx, 4), dtype=np.float64)
    # sum r*gi = dx*Su(gi*fx) - dx*Sr(gi*fy) + si*sum(gi*fval)
    for i, gi in enumerate(design):
        sb[..., i] = dx * su(gi * fx) - dx * sr(gi * fy)
        if not offset_form:
            sb[..., i] += si * _box(gi * fval, window)

    # sum r and sum r^2
    sum_r = dx * su(fx) - dx * sr(fy)
    sum_r2 = (dx * dx * (_box(c2 * fx * fx, window) - 2 * c_grid * _box(c_grid * fx * fx, window)
                         + c2 * _box(fx * fx, window))
              + dx * dx * (_box(r2 * fy * fy, window) - 2 * r_grid * _box(r_grid * fy * fy, window)
                           + r2 * _box(fy * fy, window))
              - 2 * dx * dx * (_box(cr * fx * fy, window) - c_grid * _box(r_grid * fx * fy, window)
                               - r_grid * _box(c_grid * fx * fy, window)
                               + cr * _box(fx * fy, window)))
    if not offset_form:
        sum_r = sum_r + si * _box(fval, window)
        sum_r2 = (sum_r2
                  + 2 * si * (dx * su(fx * fval) - dx * sr(fy * fval))
                  + si * si * _box(fval * fval, window))

    # ---- local strike direction from the horizontal-gradient second moment --------------
    n = _box(np.ones_like(fx), window)
    sxx = _box(fx * fx, window) / np.maximum(n, 1e-9)
    syy = _box(fy * fy, window) / np.maximum(n, 1e-9)
    sxy = _box(fx * fy, window) / np.maximum(n, 1e-9)
    tr = sxx + syy
    disc = np.sqrt(np.maximum(tr * tr / 4.0 - (sxx * syy - sxy * sxy), 0.0))
    lmaj = np.maximum(tr / 2.0 + disc, 1e-30)      # across-strike
    lmin = np.maximum(tr / 2.0 - disc, 0.0)        # along-strike
    ratio = lmin / lmaj
    ang = 0.5 * np.arctan2(2.0 * sxy, sxx - syy)   # angle of the LARGEST eigenvector
    strike_e = -np.sin(ang)                        # rotate by 90 deg -> smallest eigenvector
    strike_n = np.cos(ang)
    # A *small* minor/major ratio means one dominant gradient direction, i.e. a clean
    # straight contact whose along-strike coordinate is unidentifiable -> full penalty.
    # A ratio near 1 means an isotropic gradient (a compact source), where every
    # coordinate is identifiable -> the penalty must vanish.
    lo, hi = cfg.isotropic_fade
    fade = np.clip((hi - ratio) / max(hi - lo, 1e-9), 0.0, 1.0)

    # ---- scale each column by its *per-window* RMS for a well-conditioned solve ----------
    # A global scale would leave windows far from any anomaly badly conditioned, because the
    # normal-matrix entries there are orders of magnitude smaller than the grid average.
    n_safe = np.maximum(n, 1e-9)
    scales = np.stack([np.sqrt(np.maximum(_box(gi * gi, window) / n_safe, 1e-300))
                       for gi in design], axis=-1)          # (ny, nx, 4)
    a_sc = a / (scales[..., :, None] * scales[..., None, :])
    b_sc = sb / scales

    # ---- anisotropic regularisation: floor the along-strike constraining power -----------
    # ``strike`` is the direction of *least* horizontal-gradient variance, i.e. along the
    # contact; across-strike is perpendicular to it.  Both powers are taken from the
    # horizontal 2x2 block of the *unscaled* normal matrix, so their ratio is dimensionless
    # and independent of the units of the field.
    a_xy = a[..., :2, :2]
    across_e, across_n = -strike_n, strike_e
    power_along = (a_xy[..., 0, 0] * strike_e * strike_e
                   + 2.0 * a_xy[..., 0, 1] * strike_e * strike_n
                   + a_xy[..., 1, 1] * strike_n * strike_n)
    power_across = (a_xy[..., 0, 0] * across_e * across_e
                    + 2.0 * a_xy[..., 0, 1] * across_e * across_n
                    + a_xy[..., 1, 1] * across_n * across_n)
    lam = fade * np.maximum(cfg.min_along_strike_power * power_across - power_along, 0.0)
    # penalty lam * (theta . strike)^2 ; convert to the scaled parameterisation
    # theta_scaled_i = theta_i * scale_i
    pe = strike_e / scales[..., 0]
    pn = strike_n / scales[..., 1]
    a_sc[..., 0, 0] = a_sc[..., 0, 0] + lam * pe * pe
    a_sc[..., 1, 1] = a_sc[..., 1, 1] + lam * pn * pn
    a_sc[..., 0, 1] = a_sc[..., 0, 1] + lam * pe * pn
    a_sc[..., 1, 0] = a_sc[..., 0, 1]

    # ---- subsample to the stride lattice, then solve -----------------------------------
    half = window // 2
    rows_idx = np.asarray(rows_tops, dtype=int) - row_offset
    cols_idx = np.arange(half, nx - half, stride)
    sel = (rows_idx[:, None], cols_idx[None, :])
    sel_ix = np.ix_(rows_idx, cols_idx)

    a_sub = np.ascontiguousarray(a_sc[sel_ix])
    b_sub = np.ascontiguousarray(b_sc[sel_ix])
    eye = np.eye(4)[None, None, :, :] * 1e-10
    inv_pen = np.linalg.inv(a_sub + eye)
    # theta_scaled = theta * scale, so recover theta by dividing, not multiplying.
    sol_sc = np.einsum("...ij,...j->...i", inv_pen, b_sub)
    sol = sol_sc / scales[sel_ix]
    x0, y0, z0 = sol[..., 0], sol[..., 1], sol[..., 2]

    # ---- misfit ------------------------------------------------------------------------
    # RSS = sum r^2 - 2 s.b + s.A_unpenalised.s   (A and b are the raw window sums)
    a_raw = np.ascontiguousarray(a[sel_ix])
    b_raw = np.ascontiguousarray(sb[sel_ix])
    quad = np.einsum("...i,...ij,...j->...", sol, a_raw, sol)
    sr2 = np.maximum(sum_r2[sel_ix], 0.0)
    srn = sum_r[sel_ix]
    nn = np.maximum(n[sel_ix], 1e-9)
    rss = np.maximum(sr2 - 2.0 * np.einsum("...i,...i->...", sol, b_raw) + quad, 0.0)
    tss = np.maximum(sr2 - srn * srn / nn, 1e-300)
    rel_resid = np.sqrt(rss / tss)

    dof = np.maximum(nn - 4.0, 1.0)
    sigma2 = rss / dof
    # conditional standard error of the depth from the *penalised* covariance: the ridge is
    # part of the estimator, and the along-strike component is genuinely not identified from
    # the data, so this is a model-conditional error, not a full posterior.
    var_z = inv_pen[..., 2, 2] * (scales[sel_ix][..., 2] ** 2)
    depth_se = np.sqrt(np.maximum(var_z * sigma2, 0.0))

    offset = np.hypot(x0, y0)
    along = x0 * strike_e[sel_ix] + y0 * strike_n[sel_ix]
    cross = -x0 * strike_n[sel_ix] + y0 * strike_e[sel_ix]

    out = dict(
        row_idx=rows_idx, col_idx=cols_idx,
        centre_col=np.broadcast_to(cols_idx.astype(np.float64)[None, :], x0.shape),
        centre_row=np.broadcast_to(rows_idx.astype(np.float64)[:, None], x0.shape),
        x0=x0, y0=y0, z0=z0, rel_resid=rel_resid, depth_se=depth_se,
        offset=offset, along=along, cross=cross, ratio=ratio[sel_ix],
        n_valid=n_valid[sel_ix], last=sol[..., 3],
    )
    return out


def run_family(values: np.ndarray, valid: np.ndarray, fam: Family,
               cfg: SolveConfig) -> SolutionCloud:
    """Run one Euler family and apply the quality-control gates."""
    fld, fx, fy, fz = spectral_derivatives(values, valid, cfg.dx,
                                           upward_m=fam.upward_m,
                                           order=fam.derivative_order)
    n_valid = uniform_filter(valid.astype(np.float64), size=fam.window,
                             mode="constant", cval=0.0) * float(fam.window ** 2)
    res = solve_euler_grid(fx, fy, fz, fld, fam.si, fam.window, fam.stride, cfg, n_valid)

    ny, nx = valid.shape
    counts: dict[str, int] = {"windows": int(res["z0"].size)}

    ok = (np.isfinite(res["z0"]) & np.isfinite(res["rel_resid"]) & np.isfinite(res["depth_se"])
          & (res["n_valid"] >= cfg.min_window_valid_frac * fam.window ** 2))
    counts["rejected_incomplete_or_nonfinite"] = int((~ok).sum())
    ok &= (res["z0"] >= cfg.min_depth_m) & (res["z0"] <= cfg.max_depth_m)
    counts["after_depth_range"] = int(ok.sum())
    ok &= res["rel_resid"] <= cfg.max_rel_residual
    counts["after_rel_residual"] = int(ok.sum())
    ok &= res["depth_se"] <= np.maximum(cfg.max_depth_se_frac * np.abs(res["z0"]), 50.0)
    counts["after_depth_se"] = int(ok.sum())
    ok &= res["offset"] <= cfg.max_abs_offset_frac * 0.5 * fam.window * cfg.dx
    counts["after_offset"] = int(ok.sum())

    # suppress windows near the rectangular grid edge or near invalid data
    rows_idx, cols_idx = res["row_idx"], res["col_idx"]
    edge = np.zeros((ny, nx), bool)
    edge[:cfg.edge_suppress_cells] = True
    edge[-cfg.edge_suppress_cells:] = True
    edge[:, :cfg.edge_suppress_cells] = True
    edge[:, -cfg.edge_suppress_cells:] = True
    near_bad = uniform_filter((~valid).astype(np.float64), size=fam.window,
                              mode="constant", cval=0.0) > 0.0
    blocked = edge[np.ix_(rows_idx, cols_idx)] | near_bad[np.ix_(rows_idx, cols_idx)]
    ok &= ~blocked
    counts["after_edge_and_invalid_suppression"] = int(ok.sum())
    counts["windows_blocked_by_edge_or_invalid"] = int(blocked.sum())

    rr, cc = np.nonzero(ok)
    col = res["centre_col"][rr, cc] + res["x0"][rr, cc] / cfg.dx
    row = res["centre_row"][rr, cc] - res["y0"][rr, cc] / cfg.dx   # y north -> row decreases

    return SolutionCloud(
        family=fam.name, kind=fam.kind, si=fam.si, window=fam.window,
        col=col, row=row, depth=res["z0"][rr, cc], depth_se=res["depth_se"][rr, cc],
        rel_residual=res["rel_resid"][rr, cc], offset=res["offset"][rr, cc],
        along_strike=res["along"][rr, cc], cross_strike=res["cross"][rr, cc],
        strike_ratio=res["ratio"][rr, cc], n_obs=res["n_valid"][rr, cc],
        meta=dict(counts=counts, dx=cfg.dx,
                  derivative_order=fam.derivative_order, upward_m=fam.upward_m),
    )
