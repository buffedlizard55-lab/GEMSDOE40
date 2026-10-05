"""Euler deconvolution of potential-field data, with depth-aware solution clustering.

Implements the sliding-window Euler deconvolution of

    Reid, A. B., Allsop, J. M., Granser, H., Millett, A. J., and Somerton, I. W. (1990),
    "Magnetic interpretation in three dimensions using Euler deconvolution", Geophysics 55(1), 80-91.
    DOI 10.1190/1.1442774  (publisher: https://doi.org/10.1190/1.1442774)
    Structural-index table (ibid. Table 1; Thompson 1982 Table 1): contact / fault -> SI = 0,
    dyke / sill -> 1, horizontal cylinder -> 2, sphere / dipole -> 3.

Equation (ibid. Eq. 5, the homogeneous Euler equation for a potential field T):

    (x - x0) dT/dx + (y - y0) dT/dy + (z - z0) dT/dz = N (B - T)

with N the structural index (the negative of the homogeneity degree, SI in the geophysical
literature) and B the regional (background) field level.  Rearranged for a *measured plane* at
z = 0 and a source at depth z0 > 0 below that plane, and collecting the unknowns
[x0, y0, z0, B] linearly:

    x0*Tx + y0*Ty + z0*Tz + N*B = x*Tx + y*Ty + N*T                      (linear system, per window)

Every window is solved by ordinary least squares over the window's pixels.  The system matrix,
its normal matrix and the right-hand side are evaluated with exact box sums (``scipy.ndimage.
uniform_filter``) so the whole grid is inverted vectorially -- no per-window Python loop.

Derivative convention used here (Blakely 1995, "Potential Theory in Gravity and Magnetic
Applications", Cambridge UP, §12.5; ISBN 9780521415088): for a field known on the plane z = 0,

    Tx = F^-1[i kx F[T]],   Ty = F^-1[i ky F[T]],   Tz = F^-1[|k| F[T]]

with z positive DOWN, so that a solution at (x0, y0, z0) has z0 = depth of the contact below
the survey plane.  ``tests/test_euler.py`` verifies the sign and the depth accuracy against a
synthetic source of analytically known depth.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import distance_transform_edt, gaussian_filter, uniform_filter

# ---------------------------------------------------------------------------------------------
# 1. Grid preparation and Fourier derivatives
# ---------------------------------------------------------------------------------------------


def fill_invalid(a: np.ndarray, invalid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Nearest-neighbour fill of invalid (nodata) pixels; returns (filled, valid mask)."""
    valid = ~invalid & np.isfinite(a)
    if valid.all():
        return np.asarray(a, dtype=np.float64), valid
    idx = distance_transform_edt(~valid, return_distances=False, return_indices=True)
    filled = np.asarray(a, dtype=np.float64)[idx[0], idx[1]]
    return filled, valid


def tukey_taper(n: int, alpha: float = 0.05) -> np.ndarray:
    """Tukey (tapered cosine) window, used to suppress FFT wrap-around at the grid edges."""
    if alpha <= 0:
        return np.ones(n)
    w = np.ones(n)
    edge = max(int(alpha * n), 1)
    t = np.linspace(0.0, 1.0, edge, endpoint=False)
    ramp = 0.5 * (1.0 - np.cos(np.pi * t))
    w[:edge] = ramp
    w[-edge:] = ramp[::-1]
    return w


def fourier_gradients(grid: np.ndarray, cell_m: float = 100.0,
                      upward_continuation_m: float = 0.0,
                      lowpass_wavelength_m: float | None = None,
                      taper_alpha: float = 0.03) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (Tx, Ty, Tz) of ``grid`` on the plane z = 0, z positive down.

    Parameters
    ----------
    upward_continuation_m : float
        Continue the field upward by this height before differentiating.  This is standard
        Euler practice: it attenuates the shallowest (noise-dominated) wavelengths.
    lowpass_wavelength_m : float or None
        Gaussian low-pass cut-off wavelength.  ``None`` disables the low-pass.
    """
    g = np.asarray(grid, dtype=np.float64)
    ny, nx = g.shape
    wy = tukey_taper(ny, taper_alpha)
    wx = tukey_taper(nx, taper_alpha)
    g = g * wy[:, None] * wx[None, :]

    ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=cell_m)
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=cell_m)
    KX, KY = np.meshgrid(kx, ky)
    KR = np.hypot(KX, KY)

    G = np.fft.fft2(g)
    if lowpass_wavelength_m is not None:
        kc = 2.0 * np.pi / float(lowpass_wavelength_m)
        G = G * np.exp(-0.5 * (KR / kc) ** 2)
    if upward_continuation_m:
        G = G * np.exp(-KR * float(upward_continuation_m))

    Tx = np.real(np.fft.ifft2(1j * KX * G))
    Ty = np.real(np.fft.ifft2(1j * KY * G))
    Tz = np.real(np.fft.ifft2(KR * G))
    return Tx, Ty, Tz


# ---------------------------------------------------------------------------------------------
# 2. Sliding-window Euler deconvolution
# ---------------------------------------------------------------------------------------------


@dataclass
class EulerSolutions:
    """A dense cloud of depth-labelled Euler solutions (one per evaluated window)."""

    row: np.ndarray        # window-centre row index (grid pixels)
    col: np.ndarray        # window-centre column index
    x0: np.ndarray         # estimated source easting (m)
    y0: np.ndarray         # estimated source northing (m)
    depth_m: np.ndarray    # estimated source depth below the survey plane (m, > 0 = below)
    background: np.ndarray  # fitted regional level B
    rms: np.ndarray        # residual RMS of the window's LS fit (same units as the field)
    amp: np.ndarray        # RMS gradient amplitude inside the window (field units / m)
    offset_px: np.ndarray  # horizontal distance between window centre and (x0, y0), in pixels
    field: str             # which layer the cloud came from (e.g. "rtp")

    def __len__(self) -> int:
        return int(self.depth_m.size)


def euler_solutions(T: np.ndarray, tx: np.ndarray, ty: np.ndarray, tz: np.ndarray,
                    transform, si: float = 0.0, window_px: int = 11, stride: int = 2,
                    valid: np.ndarray | None = None, field: str = "?",
                    taper_alpha: float = 0.03) -> EulerSolutions:
    """Solve the linear Euler system in every window (window_px x window_px, stride px apart).

    ``si`` is the structural index N: 0 for a fault-like contact (Reid et al. 1990, Table 1).
    """
    T = np.asarray(T, dtype=np.float64)
    H, W = T.shape
    if valid is None:
        valid = np.ones((H, W), bool)

    # subsampling lattice (stride) so the inversion is affordable and the solutions are 200 m apart
    gy = np.arange(0, H, stride)
    gx = np.arange(0, W, stride)
    Ts = T[::stride, ::stride]
    txs, tys, tzs = tx[::stride, ::stride], ty[::stride, ::stride], tz[::stride, ::stride]
    Hs, Ws = Ts.shape

    # window size in subsampled pixels (odd, >= 3)
    w = max(int(round(window_px / stride)), 3)
    if w % 2 == 0:
        w += 1

    # global pixel coordinates of each subsampled node, in metres
    yy, xx = np.meshgrid(gy, gx, indexing="ij")
    X = transform.c + (xx + 0.5) * transform.a
    Y = transform.f + (yy + 0.5) * transform.e

    b = X * txs + Y * tys + si * Ts                      # right-hand side
    fields = {
        "TxTx": txs * txs, "TxTy": txs * tys, "TxTz": txs * tzs, "Tx": si * txs,
        "TyTy": tys * tys, "TyTz": tys * tzs, "Ty": si * tys,
        "TzTz": tzs * tzs, "Tz": si * tzs, "N1": np.full(Ts.shape, si * si, dtype=np.float32),
        "Txb": txs * b, "Tyb": tys * b, "Tzb": tzs * b, "Nb": si * b, "bb": b * b,
    }
    # float32 filtering keeps the 15 filtered fields inside the memory budget of the sandbox
    F = {k: uniform_filter(v.astype(np.float32), size=w, mode="nearest") for k, v in fields.items()}
    del fields
    bb2 = F["bb"].ravel().astype(np.float64)
    grad_energy = F["TxTx"] + F["TyTy"] + F["TzTz"]
    amp = np.sqrt(np.maximum(uniform_filter(grad_energy, size=w, mode="nearest"), 0.0))
    del grad_energy
    amp = amp.astype(np.float64)

    n_px = float(w * w)

    # window validity: (a) the window must not touch the Fourier taper zone, and (b) every pixel
    # of the window must be inside the valid domain (measured failure mode: taper-zone windows
    # produce spurious *shallow* solutions with tiny amplitudes -- see tests/test_euler.py)
    margin = int(np.ceil(taper_alpha * min(H, W))) + window_px
    edge_ok = np.zeros((Hs, Ws), bool)
    m = max(margin // stride, 1)
    edge_ok[m:-m, m:-m] = True

    if not valid.all():
        frac_s = uniform_filter(valid.astype(np.float64), size=w, mode="nearest")[::stride, ::stride]
        frac_s = frac_s[:Hs, :Ws]
        ok_win = (frac_s >= 0.999999) & edge_ok
    else:
        ok_win = edge_ok

    # regularise the normal matrix before inversion (Tikhonov on the 4 unknowns); the design
    # columns have wildly different units (nT/m vs nT), so an unregularised solve is unstable.
    scale = np.array([[max(float(np.abs(F["TxTx"]).max()), 1e-30), max(float(np.abs(F["TxTy"]).max()), 1e-30),
                       max(float(np.abs(F["TxTz"]).max()), 1e-30), max(float(np.abs(F["Tx"]).max()), 1e-30)],
                      [max(float(np.abs(F["TxTy"]).max()), 1e-30), max(float(np.abs(F["TyTy"]).max()), 1e-30),
                       max(float(np.abs(F["TyTz"]).max()), 1e-30), max(float(np.abs(F["Ty"]).max()), 1e-30)],
                      [max(float(np.abs(F["TxTz"]).max()), 1e-30), max(float(np.abs(F["TyTz"]).max()), 1e-30),
                       max(float(np.abs(F["TzTz"]).max()), 1e-30), max(float(np.abs(F["Tz"]).max()), 1e-30)],
                      [max(float(np.abs(F["Tx"]).max()), 1e-30), max(float(np.abs(F["Ty"]).max()), 1e-30),
                       max(float(np.abs(F["Tz"]).max()), 1e-30), max(float(np.abs(F["N1"]).max()), 1e-30)]])
    lam = 1e-10

    idx = np.flatnonzero(ok_win.ravel())
    sol = np.full((ok_win.size, 4), np.nan, dtype=np.float64)
    rms = np.full(ok_win.size, np.nan, dtype=np.float64)
    flat = {k: v.ravel() for k, v in F.items()}
    chunk = 150_000
    for i in range(0, idx.size, chunk):
        j = idx[i:i + chunk]
        m_ = j.size
        ATA = np.empty((m_, 4, 4), dtype=np.float64)
        g = lambda k: flat[k][j].astype(np.float64)
        ATA[:, 0, 0] = g("TxTx"); ATA[:, 0, 1] = ATA[:, 1, 0] = g("TxTy")
        ATA[:, 0, 2] = ATA[:, 2, 0] = g("TxTz"); ATA[:, 0, 3] = ATA[:, 3, 0] = g("Tx")
        ATA[:, 1, 1] = g("TyTy"); ATA[:, 1, 2] = ATA[:, 2, 1] = g("TyTz")
        ATA[:, 1, 3] = ATA[:, 3, 1] = g("Ty")
        ATA[:, 2, 2] = g("TzTz"); ATA[:, 2, 3] = ATA[:, 3, 2] = g("Tz")
        ATA[:, 3, 3] = g("N1")
        ATb = np.stack([g("Txb"), g("Tyb"), g("Tzb"), g("Nb")], axis=-1)
        ATA = ATA + lam * scale[None, :, :]
        x = np.linalg.solve(ATA, ATb[..., None])[..., 0]
        sol[j] = x
        rss = np.maximum(bb2[j] - np.einsum("ij,ij->i", x, ATb), 0.0)
        rms[j] = np.sqrt(rss / n_px)
    del flat, F

    x0 = sol[:, 0].reshape(Hs, Ws)
    y0 = sol[:, 1].reshape(Hs, Ws)
    depth = sol[:, 2].reshape(Hs, Ws)
    bg = sol[:, 3].reshape(Hs, Ws)
    rms = rms.reshape(Hs, Ws)

    off = np.hypot(y0 - Y, x0 - X) / abs(transform.a)
    keep = ok_win & np.isfinite(depth)
    return EulerSolutions(
        row=yy[keep], col=xx[keep], x0=x0[keep], y0=y0[keep], depth_m=depth[keep],
        background=bg[keep], rms=rms[keep], amp=amp[keep], offset_px=off[keep], field=field,
    )


# ---------------------------------------------------------------------------------------------
# 3. Depth-aware, tightness-weighted clustering of the solution cloud
# ---------------------------------------------------------------------------------------------


def cluster_weights(sol: EulerSolutions, max_depth_m: float = 4000.0,
                    max_offset_px: float = 5.0, cluster_radius_px: float = 2.5,
                    depth_scale_m: float = 300.0, min_neighbours: int = 6,
                    shallow_decay_m: float = 900.0, mad_scale_m: float = 180.0,
                    fit_scale_pct: float = 75.0, offset_sigma_px: float = 2.0,
                    amp_floor_frac: float = 0.25, amp_ref_pct: float = 90.0,
                    min_amp: float = 0.0, grid_shape: tuple[int, int] | None = None) -> dict:
    """Weight each solution by shallowness x depth-tightness x fit quality x lateral offset.

    The brief requires that "tight clusters of shallow, mutually-consistent solutions score
    higher than scattered or deep ones".  That is implemented as the product

        w = exp(-depth / shallow_decay)                        shallowness
          * clip(n_neighbours / min_neighbours, 0, 1)          solution density (tightness)
          * exp(-pooled_depth_std / mad_scale)                 mutual consistency of depth
          * exp(-0.5 (offset_px / offset_sigma)^2)             solution lies inside its own window
          * 1 / (1 + (rms/amp / ref)^2)                        quality of the window's LS fit

    All statistics are pooled over a monotone cell neighbourhood on a lattice of
    ``cluster_radius_px``-sized cells, so the computation is fully vectorised (no KD-tree) and the
    definition is auditable: a "cluster" is the set of solutions falling in the same or an
    immediately adjacent cell (8-neighbourhood), and the depth spread is the pooled standard
    deviation of those solutions' estimated depths.
    """
    n = len(sol)
    keys = dict(weight=np.zeros(0), n_neighbours=np.zeros(0), cluster_depth_mad=np.zeros(0),
                fit_weight=np.zeros(0), depth_weight=np.zeros(0), offset_weight=np.zeros(0),
                keep=np.zeros(0, bool), amp=np.zeros(0), rel_fit=np.zeros(0),
                amp_gate=np.zeros(0, bool), pooled_depth_mean=np.zeros(0),
                pooled_depth_std=np.zeros(0), cell_depth_std=np.zeros(0))
    if n == 0:
        return keys

    H, W = grid_shape if grid_shape is not None else (int(sol.row.max()) + 2, int(sol.col.max()) + 2)
    cell = max(int(np.ceil(cluster_radius_px)), 1)
    ny, nx = int(H // cell) + 2, int(W // cell) + 2
    ci = np.clip((sol.row // cell).astype(np.int64), 0, ny - 1)
    cj = np.clip((sol.col // cell).astype(np.int64), 0, nx - 1)
    flat_cell = ci * nx + cj
    depth = np.nan_to_num(sol.depth_m, nan=0.0, posinf=0.0, neginf=0.0)

    cnt = np.bincount(flat_cell, minlength=ny * nx).astype(np.float64)
    s1 = np.bincount(flat_cell, weights=depth, minlength=ny * nx)
    s2 = np.bincount(flat_cell, weights=depth * depth, minlength=ny * nx)

    def _neigh(a):
        """Sum a per-cell field over the 3x3 cell neighbourhood."""
        g = a.reshape(ny, nx)
        out = np.zeros_like(g)
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                out += np.roll(np.roll(g, dy, axis=0), dx, axis=1)
        return out

    n_nb = _neigh(cnt)
    s1_nb = _neigh(s1)
    s2_nb = _neigh(s2)
    safe_n = np.maximum(n_nb, 1.0)
    pooled_mean = s1_nb / safe_n
    pooled_var = np.maximum(s2_nb / safe_n - pooled_mean ** 2, 0.0)
    pooled_std = np.sqrt(pooled_var)

    # per-cell (not neighbourhood-pooled) spread: the depth scatter inside one cell
    safe_c = np.maximum(cnt, 1.0)
    cm = s1 / safe_c
    cvar = np.maximum(s2 / safe_c - cm ** 2, 0.0)
    cell_std = np.sqrt(cvar)

    n_i = n_nb.ravel()[flat_cell]
    std_i = pooled_std.ravel()[flat_cell]
    mean_i = pooled_mean.ravel()[flat_cell]
    cstd_i = cell_std.ravel()[flat_cell]

    density_w = np.clip(n_i / max(min_neighbours, 1), 0.0, 1.0)
    depth_w = np.exp(-np.clip(depth, 0, None) / shallow_decay_m)
    mad_w = np.exp(-np.clip(std_i, 0, None) / mad_scale_m)
    offset_w = np.exp(-0.5 * (sol.offset_px / max(offset_sigma_px, 1e-6)) ** 2)

    # fit quality relative to the field variation in the same window: an absolute residual would
    # reward flat, signal-free areas (measured failure mode; see tests/test_euler.py)
    amp = np.nan_to_num(sol.amp, nan=0.0, posinf=0.0)
    rel = np.nan_to_num(sol.rms, nan=np.inf, posinf=np.inf) / np.maximum(amp, 1e-30)
    r = rel[np.isfinite(rel)]
    ref = np.percentile(r, fit_scale_pct) if r.size else 1.0
    fit_w = 1.0 / (1.0 + (rel / max(ref, 1e-12)) ** 2)

    amp_ref = float(np.percentile(amp, amp_ref_pct)) if amp.size else 0.0
    amp_gate = amp >= max(amp_floor_frac * amp_ref, min_amp)

    keep = ((depth > 0.0) & (depth <= max_depth_m) &
            (sol.offset_px <= max_offset_px) & np.isfinite(sol.rms) & amp_gate)
    w = density_w * depth_w * mad_w * offset_w * fit_w
    w = np.where(keep, w, 0.0)
    return dict(weight=w, n_neighbours=n_i, cluster_depth_mad=cstd_i, fit_weight=fit_w,
                depth_weight=depth_w, offset_weight=offset_w, keep=keep, amp=amp,
                rel_fit=rel, amp_gate=amp_gate, pooled_depth_mean=mean_i,
                pooled_depth_std=std_i, cell_depth_std=cstd_i)


def kde_field(shape: tuple[int, int], row: np.ndarray, col: np.ndarray, weight: np.ndarray,
              sigma_px: float = 1.5) -> np.ndarray:
    """Kernel-density estimate of weighted solution density per pixel.

    Weighted solutions are accumulated on the pixel lattice and convolved with an isotropic
    Gaussian of width ``sigma_px`` (in pixels; 1 px = 100 m).  ``sigma_px`` is chosen at or below
    the metric's own 3 px triangular kernel so the emission is not blurred beyond the scale at
    which the scorer can still pay for it.
    """
    acc = np.zeros(shape, dtype=np.float64)
    np.add.at(acc, (row.astype(np.int64), col.astype(np.int64)), weight)
    if sigma_px > 0:
        acc = gaussian_filter(acc, sigma=sigma_px, mode="constant")
    return acc
