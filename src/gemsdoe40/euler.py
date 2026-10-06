"""3-D Euler deconvolution of potential-field grids (Reid et al. 1990).

Reference (must-cite, freely mirrored by the author):
    Reid, A.B., Allsop, J.M., Granser, H., Millett, A.J. and Somerton, I.W., 1990.
    Magnetic interpretation in three dimensions using Euler deconvolution.
    Geophysics, 55(1), 80–91.  https://doi.org/10.1190/1.1442774
    Author PDF: https://www.reid-geophys.co.uk/wp-content/uploads/2017/11/Reid-et-al-1990.pdf

Euler's homogeneity relation on a potential field T observed at (x, y, z):

    (x - x0) ∂T/∂x + (y - y0) ∂T/∂y + (z - z0) ∂T/∂z = N (B - T)

N is the structural index.  For a magnetic contact of great depth extent Reid
et al. 1990 (Appendix) show N = 0; for a gravity fault/contact Reid (2003)
likewise uses N = 0.  Those are the indices this module uses for a fault-like
contact.  The fourth unknown B (background) is solved jointly.

This is a *source locator*, not an edge detector: the output of each window is
a depth-labelled point (x0, y0, z0), not a gradient-magnitude pixel.  Tight
clusters of shallow, mutually-consistent points are the signature of a real
near-surface contact; scattered or deep points are the signature of noise.

Implementation notes (original to GEMSDOE40, not a port of any prior GEMSDOE
Euler arm):

  * Two independent fields are deconvolved: reduced-to-pole magnetics (band
    ``rtp``) and isostatic gravity (band ``iso_grav_anom``).  GEMSDOE28's
    H38-1 arm deconvolved TMI only and then *added 200 dots onto an H19-5
    dotted base*; this module never touches that base.
  * Horizontal derivatives are central differences in field-units per metre.
    The vertical derivative is a single-grid FFT |k|T (Blakely 1995, eq. 12-8),
    not a tiled Tukey scheme.
  * Windows are accepted only where the analytic-signal amplitude is locally
    elevated (located Euler; Salem & Ravat 2003 style gating) AND the 4×4
    normal matrix is well-conditioned AND the solved source sits inside the
    window AND relative depth uncertainty is small.
  * Two window sizes (8 px = 800 m, 12 px = 1.2 km) are run so a single depth
    family cannot dominate.

No sklearn.  No DBSCAN.  Clustering is a later KDE stage (``kde.py``).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.fft import irfft2, rfft2

PIXEL_M = 100.0


@dataclass
class EulerCloud:
    """Depth-labelled Euler solutions from one field / window / SI setting."""

    field_name: str
    structural_index: float
    window_px: int
    stride_px: int
    row: np.ndarray
    col: np.ndarray
    depth_m: np.ndarray
    rel_se: np.ndarray
    cond: np.ndarray
    analytic: np.ndarray
    stats: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return int(self.row.size)


def _fill_nearest(arr: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Nearest-valid fill so FFT / differences have context. Filled cells stay invalid."""
    from scipy.ndimage import distance_transform_edt

    if valid.all():
        return np.asarray(arr, dtype=np.float64)
    idx = distance_transform_edt(~valid, return_distances=False, return_indices=True)
    return np.asarray(arr, dtype=np.float64)[tuple(idx)]


def derivatives(field: np.ndarray, valid: np.ndarray, pixel_m: float = PIXEL_M):
    """Return (gx_east, gy_north, gz_down, analytic_signal, deriv_valid) in field-units / pixel.

    Horizontal: central differences, valid only when both neighbours are valid.
    Vertical: FFT |k| T on a nearest-filled, mean-removed grid (Blakely 1995),
    scaled by pixel_m so it is per-pixel like gx/gy.
    Analytic signal A = sqrt(gx² + gy² + gz²) (Roest, Verhoef & Pilkington 1992).
    """
    f = np.asarray(field, dtype=np.float64)
    v = np.asarray(valid, dtype=bool)
    h, w = f.shape
    filled = _fill_nearest(np.where(v, f, 0.0), v)

    # Derivatives in field-units *per pixel* so the Euler unknowns (x0, y0, z0)
    # stay in pixel units.  Depth in metres is z0_px * pixel_m after the solve.
    gx = np.zeros((h, w), dtype=np.float64)
    gy = np.zeros((h, w), dtype=np.float64)
    hvalid = np.zeros((h, w), dtype=bool)
    if w >= 3 and h >= 3:
        gx[:, 1:-1] = (filled[:, 2:] - filled[:, :-2]) * 0.5
        gy[1:-1, :] = -(filled[2:, :] - filled[:-2, :]) * 0.5
        xok = v[:, 1:-1] & v[:, :-2] & v[:, 2:]
        yok = v[1:-1, :] & v[:-2, :] & v[2:, :]
        hvalid[1:-1, 1:-1] = xok[1:-1, :] & yok[:, 1:-1]

    # Whole-grid vertical derivative.  |k|T is in field-units / metre; multiply
    # by pixel_m to convert to field-units / pixel (Blakely 1995, eq. 12-8).
    mean = float(filled[v].mean()) if v.any() else 0.0
    work = filled - mean
    ky = 2.0 * np.pi * np.fft.fftfreq(h, d=pixel_m)
    kx = 2.0 * np.pi * np.fft.rfftfreq(w, d=pixel_m)
    k = np.hypot(ky[:, None], kx[None, :])
    spec = rfft2(work)
    gz = irfft2(spec * k, s=(h, w)).real * pixel_m
    dvalid = hvalid & v
    gx = np.where(dvalid, gx, 0.0)
    gy = np.where(dvalid, gy, 0.0)
    gz = np.where(dvalid, gz, 0.0)
    analytic = np.sqrt(gx * gx + gy * gy + gz * gz)
    analytic = np.where(dvalid, analytic, 0.0)
    return gx, gy, gz, analytic, dvalid


def _box_sum(values: np.ndarray, window: int, stride: int) -> np.ndarray:
    """Summed-area-table box sums, sampled every ``stride`` pixels at the top-left."""
    a = np.asarray(values, dtype=np.float64)
    h, w = a.shape
    if window > h or window > w:
        return np.empty((0, 0), dtype=np.float64)
    sat = np.zeros((h + 1, w + 1), dtype=np.float64)
    sat[1:, 1:] = np.cumsum(np.cumsum(a, axis=0), axis=1)
    s = (
        sat[window:, window:]
        - sat[:-window, window:]
        - sat[window:, :-window]
        + sat[:-window, :-window]
    )
    return s[::stride, ::stride].copy()


def deconvolve(
    field: np.ndarray,
    valid: np.ndarray,
    *,
    field_name: str,
    structural_index: float = 0.0,
    window_px: int = 10,
    stride_px: int = 4,
    pixel_m: float = PIXEL_M,
    analytic_percentile: float = 70.0,
    max_condition: float = 5e4,
    max_rel_se: float = 0.20,
    min_depth_m: float = 50.0,
    max_depth_m: float = 2500.0,
    source_pad_px: float = 1.5,
) -> EulerCloud:
    """Moving-window least-squares Euler deconvolution (Reid et al. 1990).

    For SI = 0 the fourth unknown is the background offset A (Reid's B).
    For SI > 0 the SI·T term is moved onto the right-hand side.

    Coordinates: column increases east, row increases south.  Depth is positive
    downward in metres.  Only fully-valid windows whose mean analytic-signal
    amplitude exceeds ``analytic_percentile`` of valid cells are solved.
    """
    if structural_index < 0:
        raise ValueError("structural index must be >= 0")
    if window_px < 4 or stride_px < 1:
        raise ValueError("window_px >= 4 and stride_px >= 1 required")

    gx, gy, gz, analytic, dvalid = derivatives(field, valid, pixel_m)
    t = np.where(dvalid, np.asarray(field, dtype=np.float64), 0.0)
    h, w = t.shape
    empty = EulerCloud(
        field_name, structural_index, window_px, stride_px,
        np.empty(0), np.empty(0), np.empty(0), np.empty(0), np.empty(0), np.empty(0),
        stats={"accepted_solutions": 0, "reason": "grid smaller than window"},
    )
    if window_px > h or window_px > w:
        return empty

    n_obs = window_px * window_px
    # Gating: a window must be fully valid AND sit on an elevated analytic signal.
    win_valid_n = _box_sum(dvalid.astype(np.float64), window_px, stride_px)
    fully = win_valid_n == float(n_obs)
    win_as = _box_sum(analytic, window_px, stride_px) / float(n_obs)
    if dvalid.any():
        thr = float(np.percentile(analytic[dvalid], analytic_percentile))
    else:
        thr = np.inf
    gated = fully & (win_as >= thr)
    nr, nc = fully.shape
    stats_base = {
        "field_name": field_name,
        "structural_index": float(structural_index),
        "window_px": int(window_px),
        "stride_px": int(stride_px),
        "window_count": int(nr * nc),
        "fully_valid_windows": int(fully.sum()),
        "analytic_percentile": float(analytic_percentile),
        "analytic_threshold": float(thr) if np.isfinite(thr) else None,
        "gated_windows": int(gated.sum()),
    }
    if not gated.any():
        empty.stats = {**stats_base, "accepted_solutions": 0}
        return empty

    # Normal matrix columns: [gx, gy, gz, 1]
    moments = np.zeros((nr, nc, 4, 4), dtype=np.float64)
    cols_g = (gx, gy, gz)
    for i in range(3):
        for j in range(i, 3):
            block = _box_sum(cols_g[i] * cols_g[j], window_px, stride_px)
            moments[:, :, i, j] = block
            moments[:, :, j, i] = block
    for i, g in enumerate(cols_g):
        s = _box_sum(g, window_px, stride_px)
        moments[:, :, i, 3] = s
        moments[:, :, 3, i] = s
    moments[:, :, 3, 3] = float(n_obs)

    # Observation coordinates in pixel units, origin at grid centre so the
    # fourth column and the (x, y) unknowns stay similarly scaled.
    x_east = np.arange(w, dtype=np.float64) - 0.5 * (w - 1)
    y_north = 0.5 * (h - 1) - np.arange(h, dtype=np.float64)
    rhs_field = gx * x_east[None, :] + gy * y_north[:, None]
    if structural_index > 0:
        rhs_field = rhs_field + structural_index * t

    rhs = np.empty((nr, nc, 4), dtype=np.float64)
    rhs[:, :, 0] = _box_sum(gx * rhs_field, window_px, stride_px)
    rhs[:, :, 1] = _box_sum(gy * rhs_field, window_px, stride_px)
    rhs[:, :, 2] = _box_sum(gz * rhs_field, window_px, stride_px)
    rhs[:, :, 3] = _box_sum(rhs_field, window_px, stride_px)
    y2 = _box_sum(rhs_field * rhs_field, window_px, stride_px)

    # Column-normalise, then solve only gated windows.
    colnorm = np.sqrt(np.maximum(np.diagonal(moments, axis1=2, axis2=3), 1e-30))
    scaled = moments / (colnorm[:, :, :, None] * colnorm[:, :, None, :])
    # Eigenvalues of the 4×4 for a cheap condition-number gate.
    evals = np.linalg.eigvalsh(scaled)
    cond = np.full((nr, nc), np.inf)
    ok_e = (evals[:, :, 0] > 1e-12) & np.isfinite(evals[:, :, -1])
    cond[ok_e] = np.sqrt(evals[:, :, -1][ok_e] / evals[:, :, 0][ok_e])
    well = gated & ok_e & (cond <= max_condition)
    flat = np.flatnonzero(well.ravel())
    stats_base["condition_pass"] = int(flat.size)
    if flat.size == 0:
        empty.stats = {**stats_base, "accepted_solutions": 0}
        return empty

    nrm = colnorm.reshape(-1, 4)[flat]
    A = scaled.reshape(-1, 4, 4)[flat]
    b = rhs.reshape(-1, 4)[flat] / nrm
    try:
        beta_s = np.linalg.solve(A, b[..., None])[..., 0]
    except np.linalg.LinAlgError:
        empty.stats = {**stats_base, "accepted_solutions": 0, "reason": "linalg"}
        return empty
    beta = beta_s / nrm  # unscale → (x0_east, y0_north, z0_px, offset)
    rss = np.maximum(y2.ravel()[flat] - np.einsum("ij,ij->i", beta, rhs.reshape(-1, 4)[flat]), 0.0)
    var = rss / max(n_obs - 4, 1)
    inv = np.linalg.inv(A)
    z_var = var * inv[:, 2, 2] / np.square(nrm[:, 2])
    depth_se_px = np.sqrt(np.maximum(z_var, 0.0))
    rel_se = depth_se_px / np.maximum(np.abs(beta[:, 2]), 1e-12)
    depth_m = beta[:, 2] * pixel_m

    rows = np.arange(0, h - window_px + 1, stride_px, dtype=np.int32)
    cols = np.arange(0, w - window_px + 1, stride_px, dtype=np.int32)
    wr, wc = np.unravel_index(flat, (nr, nc))
    top_r = rows[wr]
    top_c = cols[wc]
    src_c = beta[:, 0] + 0.5 * (w - 1)
    src_r = 0.5 * (h - 1) - beta[:, 1]

    depth_ok = np.isfinite(depth_m) & (depth_m >= min_depth_m) & (depth_m <= max_depth_m)
    se_ok = np.isfinite(rel_se) & (rel_se <= max_rel_se)
    in_win = (
        (src_c >= top_c - source_pad_px)
        & (src_c <= top_c + window_px - 1 + source_pad_px)
        & (src_r >= top_r - source_pad_px)
        & (src_r <= top_r + window_px - 1 + source_pad_px)
    )
    keep = depth_ok & se_ok & in_win
    sel = np.flatnonzero(keep)

    # Analytic-signal amplitude at the solved (rounded) location, for weighting.
    rr = np.clip(np.rint(src_r[sel]).astype(int), 0, h - 1)
    cc = np.clip(np.rint(src_c[sel]).astype(int), 0, w - 1)
    as_at = analytic[rr, cc]

    stats = {
        **stats_base,
        "depth_pass": int(depth_ok.sum()),
        "rel_se_pass": int((depth_ok & se_ok).sum()),
        "accepted_solutions": int(sel.size),
        "depth_median_m": float(np.median(depth_m[sel])) if sel.size else None,
        "depth_p90_m": float(np.percentile(depth_m[sel], 90)) if sel.size else None,
    }
    return EulerCloud(
        field_name=field_name,
        structural_index=float(structural_index),
        window_px=int(window_px),
        stride_px=int(stride_px),
        row=src_r[sel].astype(np.float64),
        col=src_c[sel].astype(np.float64),
        depth_m=depth_m[sel].astype(np.float64),
        rel_se=rel_se[sel].astype(np.float64),
        cond=cond.ravel()[flat][sel].astype(np.float64),
        analytic=as_at.astype(np.float64),
        stats=stats,
    )


def merge_clouds(clouds: list[EulerCloud]) -> EulerCloud:
    """Concatenate several Euler clouds (different fields / windows) into one."""
    if not clouds:
        raise ValueError("no clouds")
    if len(clouds) == 1:
        return clouds[0]
    return EulerCloud(
        field_name="+".join(sorted({c.field_name for c in clouds})),
        structural_index=float(np.mean([c.structural_index for c in clouds])),
        window_px=-1,
        stride_px=-1,
        row=np.concatenate([c.row for c in clouds]),
        col=np.concatenate([c.col for c in clouds]),
        depth_m=np.concatenate([c.depth_m for c in clouds]),
        rel_se=np.concatenate([c.rel_se for c in clouds]),
        cond=np.concatenate([c.cond for c in clouds]),
        analytic=np.concatenate([c.analytic for c in clouds]),
        stats={"n_clouds": len(clouds), "accepted_solutions": int(sum(len(c) for c in clouds)),
               "per_cloud": [c.stats for c in clouds]},
    )
