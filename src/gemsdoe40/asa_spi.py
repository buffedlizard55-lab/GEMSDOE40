"""H8: analytic-signal / local-wavenumber source imaging (SPI family).

Physics
-------
For a harmonic potential field T the analytic-signal amplitude

    A = sqrt(Tx^2 + Ty^2 + Tz^2)

peaks directly over magnetization/density contrast boundaries independently
of the magnetization direction (Nabighian 1972; Roest, Verhoef & Pilkington
1992).  Where T has a contact-like (structural index 0) singularity at
horizontal position (x0, y0) and depth z0, A(ρ) = C (ρ² + z0²)^(−1/2) in
cross-strike distance ρ, so at the peak the Hessian eigenvalue of A is
λ_min = −C/z0³ and A = C/z0, giving the contact-calibrated local-wavenumber
depth (SPI family; Thurston & Smith 1997; Blakely 1995, ch. 3)

    z0 = sqrt(A / |λ_min|).

A source whose analytic signal decays as r^(−(n+1)) returns z_true/√(n+1):
exact for contacts (n = 0), z_true/√2 for pole-like sources (n = 1); the
synthetic suite asserts both.  The Hessian is taken by central differences
of the smooth amplitude field A; Fourier second derivatives of T are not
used (the preregistration amendment records why).

Unlike Euler deconvolution (a windowed least-squares solve), this is a
peak/derivative-ratio estimator with different failure modes, so agreement
between the two is independent evidence.  The analytic signal appears
elsewhere in this repository ONLY as a window quality gate
(``gemsdoe40.euler``); here it is the estimator itself.

The gravity family applies the same pipeline to the first vertical
derivative of the isostatic gravity anomaly as a local top-edge
approximation (Reid & Thurston 2014), exactly as the frozen H4 Euler
experiment did.

Nothing in this module reads labels, proxy truth, prior predictions or
scores.  Catalogue handling (exact-pixel mask + flank ramp) is applied by
the caller from the publicly supplied known-fault raster.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy import fft as sp_fft
from scipy import ndimage
from scipy.spatial import cKDTree

RESOLUTION_M = 100.0
PAD_CELLS = 128
EPS = 1e-12
# Supergaussian taper on the wavenumber spectrum (preregistration item-3 amendment).
# Untapered Fourier derivatives of contact-like fields carry reflect-pad seam
# ringing and aliased grid-scale oscillation that corrupt the amplitude
# curvature used for depth.  The taper kills the highest wavenumbers; the
# measured trade-off is that the shallowest ridges (z0 ~ 3 cells) are biased
# ~25% deep while z0 >= 8 cells is recovered within ~5%.
TAPER_KC_FRAC = 0.45
TAPER_POWER = 4
# Gaussian smoothing of the amplitude field before the curvature depth step:
# stabilizes the Hessian against noise while preserving ridges >= ~5 cells wide
# (measured on the registered synthetic suite).
AMPLITUDE_SMOOTH_PX = 1.5


@dataclass
class PeakCloud:
    """Depth-labelled analytic-signal peaks, coordinates in grid pixels, depth in metres."""

    family: str
    row: np.ndarray
    col: np.ndarray
    depth_m: np.ndarray
    amplitude: np.ndarray
    persistent: np.ndarray
    corroborated: np.ndarray
    weight: np.ndarray
    stats: dict[str, Any] = field(default_factory=dict)

    def __len__(self) -> int:
        return int(self.row.size)


def _finite_mask(array: np.ndarray) -> np.ndarray:
    arr = np.asarray(array)
    return np.isfinite(arr) & (np.abs(arr) < 1e30)


def fill_nearest(array: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Nearest-fill invalid cells strictly for derivative computation (float32)."""
    arr = np.asarray(array, dtype=np.float32)
    valid = _finite_mask(arr)
    if not valid.any():
        raise ValueError("cannot fill a band with no finite values")
    if valid.all():
        return arr.copy(), valid
    idx = ndimage.distance_transform_edt(~valid, return_distances=False, return_indices=True)
    filled = arr.copy()
    filled[~valid] = arr[tuple(idx)][~valid]
    del idx
    return filled, valid


def fourier_derivatives(
    array: np.ndarray,
    *,
    pad_cells: int = PAD_CELLS,
    resolution_m: float = RESOLUTION_M,
    taper_kc_frac: float = TAPER_KC_FRAC,
    taper_power: int = TAPER_POWER,
) -> dict[str, np.ndarray]:
    """Tapered FFT derivatives of a harmonic field on a padded reflect grid.

    Grid convention: rows increase southward (the sample transform has a
    -100 m north step), so the north derivative is -d/drow.  Returns
    Tx (east), Ty (north), Tz (downward positive) in field-units per metre,
    cropped back to the input shape.  All spectra are multiplied by the
    supergaussian taper documented at module scope.
    """
    filled, valid = fill_nearest(array)
    mean = float(filled[valid].mean())
    padded = np.pad((filled - np.float32(mean)).astype(np.float32),
                    ((pad_cells, pad_cells), (pad_cells, pad_cells)), mode="reflect")
    del filled
    ny, nx = padded.shape
    fy = (2.0 * np.pi * np.fft.fftfreq(ny, d=resolution_m)).astype(np.float32)[:, None]
    fx = (2.0 * np.pi * np.fft.rfftfreq(nx, d=resolution_m)).astype(np.float32)[None, :]
    k = np.sqrt(fy * fy + fx * fx).astype(np.float32)
    k_nyquist = np.pi / resolution_m
    taper = np.exp(-(k / np.float32(taper_kc_frac * k_nyquist)) ** taper_power).astype(np.float32)
    spec = sp_fft.rfft2(padded, workers=2)
    del padded
    crop = (slice(pad_cells, -pad_cells), slice(pad_cells, -pad_cells))

    def back(arr_spec: np.ndarray) -> np.ndarray:
        out = sp_fft.irfft2(arr_spec, s=(ny, nx), workers=2)[crop]
        return np.ascontiguousarray(out, dtype=np.float32)

    out = {
        "tx": back((1j * fx) * spec * taper),
        "ty": back((-1j * fy) * spec * taper),
        "tz": back(k * spec * taper),
    }
    del spec, taper, k, fx, fy
    return out


def vertical_derivative(
    array: np.ndarray,
    *,
    pad_cells: int = PAD_CELLS,
    resolution_m: float = RESOLUTION_M,
    taper_kc_frac: float = TAPER_KC_FRAC,
    taper_power: int = TAPER_POWER,
) -> np.ndarray:
    """Downward-positive tapered +|k| derivative (used to form Gz for gravity)."""
    filled, valid = fill_nearest(array)
    mean = float(filled[valid].mean())
    padded = np.pad((filled - np.float32(mean)).astype(np.float32),
                    ((pad_cells, pad_cells), (pad_cells, pad_cells)), mode="reflect")
    del filled
    ny, nx = padded.shape
    fy = (2.0 * np.pi * np.fft.fftfreq(ny, d=resolution_m)).astype(np.float32)[:, None]
    fx = (2.0 * np.pi * np.fft.rfftfreq(nx, d=resolution_m)).astype(np.float32)[None, :]
    k = np.sqrt(fy * fy + fx * fx).astype(np.float32)
    k_nyquist = np.pi / resolution_m
    taper = np.exp(-(k / np.float32(taper_kc_frac * k_nyquist)) ** taper_power).astype(np.float32)
    spec = sp_fft.rfft2(padded, workers=2)
    del padded
    out = sp_fft.irfft2(k * taper * spec, s=(ny, nx), workers=2)[pad_cells:-pad_cells, pad_cells:-pad_cells]
    del spec, k, taper, fx, fy
    return np.ascontiguousarray(out, dtype=np.float32)


def upward_continue(
    array: np.ndarray,
    height_m: float,
    *,
    pad_cells: int = PAD_CELLS,
    resolution_m: float = RESOLUTION_M,
) -> np.ndarray:
    """Continue a harmonic field upward by exp(-|k|h)."""
    if height_m <= 0:
        return np.asarray(array, dtype=np.float32)
    filled, valid = fill_nearest(array)
    mean = float(filled[valid].mean())
    padded = np.pad((filled - np.float32(mean)).astype(np.float32),
                    ((pad_cells, pad_cells), (pad_cells, pad_cells)), mode="reflect")
    del filled
    ny, nx = padded.shape
    fy = (2.0 * np.pi * np.fft.fftfreq(ny, d=resolution_m)).astype(np.float32)[:, None]
    fx = (2.0 * np.pi * np.fft.rfftfreq(nx, d=resolution_m)).astype(np.float32)[None, :]
    k = np.sqrt(fy * fy + fx * fx).astype(np.float32)
    spec = sp_fft.rfft2(padded, workers=2)
    del padded
    spec *= np.exp(-k * np.float32(height_m)).astype(np.float32)
    out = sp_fft.irfft2(spec, s=(ny, nx), workers=2)[pad_cells:-pad_cells, pad_cells:-pad_cells]
    return np.ascontiguousarray(out, dtype=np.float32)


def hessian_min_eigenvalue(a: np.ndarray, resolution_m: float = RESOLUTION_M) -> np.ndarray:
    """Smaller eigenvalue of the Hessian of A (per m^2), central differences.

    At an analytic-signal peak/ridge this is negative; its magnitude encodes
    the source depth for the contact-calibrated estimator.
    """
    ny, nx = a.shape
    lam = np.zeros((ny, nx), dtype=np.float64)
    if ny < 3 or nx < 3:
        return lam
    inv_h2 = 1.0 / (resolution_m * resolution_m)
    axx = (a[:, 2:] - 2.0 * a[:, 1:-1] + a[:, :-2]) * inv_h2
    ayy = (a[2:, :] - 2.0 * a[1:-1, :] + a[:-2, :]) * inv_h2
    axy = (a[2:, 2:] - a[2:, :-2] - a[:-2, 2:] + a[:-2, :-2]) * (0.25 * inv_h2)
    mean_c = (axx[1:-1, :] + ayy[:, 1:-1]) * 0.5
    diff_c = (axx[1:-1, :] - ayy[:, 1:-1]) * 0.5
    disc = np.sqrt(diff_c * diff_c + axy * axy)
    lam[1:-1, 1:-1] = mean_c - disc
    return lam


def analytic_signal_depth(
    array: np.ndarray,
    *,
    pad_cells: int = PAD_CELLS,
    resolution_m: float = RESOLUTION_M,
) -> dict[str, np.ndarray]:
    """Return A and the contact-calibrated depth estimate z = sqrt(A/|lambda_min|).

    The amplitude A peaks over source edges (Nabighian 1972); the Hessian of
    the smooth amplitude field supplies the local-wavenumber depth without
    Fourier second derivatives of T (see module docstring and the
    preregistration amendment).
    """
    d = fourier_derivatives(array, pad_cells=pad_cells, resolution_m=resolution_m)
    a = np.sqrt(d["tx"] ** 2 + d["ty"] ** 2 + d["tz"] ** 2)
    a_smooth = ndimage.gaussian_filter(a, AMPLITUDE_SMOOTH_PX, truncate=3.0)
    lam = hessian_min_eigenvalue(a_smooth, resolution_m)
    z = np.full(a.shape, np.inf, dtype=np.float64)
    ok = lam < 0
    z[ok] = np.sqrt(a_smooth[ok] / np.abs(lam[ok]))
    return {"a": a, "a_smooth": a_smooth, "lambda_min": lam, "z": z, **d}


def far_from_invalid(valid: np.ndarray, margin_m: float, resolution_m: float = RESOLUTION_M) -> np.ndarray:
    """Cells at least margin_m from any invalid cell or the rectangular edge."""
    dist_px = ndimage.distance_transform_edt(valid)
    edge = np.ones_like(valid)
    edge[0, :] = edge[-1, :] = edge[:, 0] = edge[:, -1] = False
    edge_dist = ndimage.distance_transform_edt(edge)
    d = np.minimum(dist_px, edge_dist)
    return d * resolution_m >= margin_m


def local_maxima(a: np.ndarray, where: np.ndarray, threshold: float) -> np.ndarray:
    """3x3 non-maximum suppression mask."""
    dil = ndimage.maximum_filter(a, size=3)
    return where & (a >= threshold) & (a == dil)


def select_peaks(
    a: np.ndarray,
    z: np.ndarray,
    valid: np.ndarray,
    *,
    margin_m: float,
    quantile: float = 0.95,
    depth_bounds_m: tuple[float, float] = (100.0, 5000.0),
    median_window: int = 5,
    median_tolerance_factor: float = 0.35,
    median_tolerance_floor_m: float = 250.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, int]]:
    """NMS peaks with amplitude, depth-QC and local median consistency gates."""
    where = valid & far_from_invalid(valid, margin_m)
    stats: dict[str, int] = {}
    threshold = float(np.quantile(a[where], quantile)) if where.any() else np.inf
    peaks = local_maxima(a, where, threshold)
    stats["peaks_nms_threshold"] = int(peaks.sum())
    row, col = np.nonzero(peaks)
    zsel = z[row, col]
    keep = (zsel >= depth_bounds_m[0]) & (zsel <= depth_bounds_m[1])
    stats["peaks_depth_qc"] = int(keep.sum())
    med = ndimage.median_filter(np.where(where, z, 0.0), size=median_window)
    tol = np.maximum(median_tolerance_floor_m, median_tolerance_factor * np.abs(zsel))
    keep &= np.abs(zsel - med[row, col]) <= tol
    stats["peaks_median_consistent"] = int(keep.sum())
    return row[keep], col[keep], zsel[keep], stats


def match_persistence(
    base_row: np.ndarray,
    base_col: np.ndarray,
    base_z: np.ndarray,
    cont_row: np.ndarray,
    cont_col: np.ndarray,
    cont_z_raw: np.ndarray,
    height_m: float,
    *,
    xy_tol_px: float = 3.0,
    depth_tol_m: float = 400.0,
    depth_rel: float = 0.4,
) -> np.ndarray:
    """True where a continued-level peak corroborates the base peak.

    The continued observation plane is farther from the source, so its raw
    contact-depth estimate exceeds the base estimate by ~height_m.
    """
    persistent = np.zeros(base_row.size, dtype=bool)
    if base_row.size == 0 or cont_row.size == 0:
        return persistent
    cont_z = cont_z_raw - height_m
    tree = cKDTree(np.column_stack([cont_row, cont_col]))
    dist, idx = tree.query(np.column_stack([base_row, base_col]), k=1, distance_upper_bound=xy_tol_px)
    found = np.isfinite(dist)
    if not found.any():
        return persistent
    tol = np.maximum(depth_tol_m, depth_rel * np.abs(base_z))
    cand = np.where(found, idx, 0)
    z_match = np.abs(cont_z[cand] - base_z) <= tol
    persistent = found & z_match
    return persistent


def cloud_weights(
    row: np.ndarray,
    col: np.ndarray,
    z: np.ndarray,
    persistent: np.ndarray,
    corroborated: np.ndarray,
    *,
    depth_scale_m: float = 1500.0,
    neighbor_px: float = 3.0,
    mad_floor_m: float = 250.0,
    mad_rel: float = 0.35,
    corroborated_boost: float = 2.5,
    non_persistent_factor: float = 0.35,
) -> np.ndarray:
    """Registered H8 weight: shallow x persistent x tight x depth-consistent x corroborated."""
    n = row.size
    if n == 0:
        return np.empty(0, dtype=np.float64)
    shallow = np.exp(-np.clip(z, 0.0, None) / depth_scale_m)
    persist = np.where(persistent, 1.0, non_persistent_factor)
    xy = np.column_stack([row, col]).astype(np.float64)
    tree = cKDTree(xy)
    pairs = tree.query_pairs(r=neighbor_px, output_type="ndarray")
    counts = np.zeros(n, dtype=np.float64)
    scatter = np.zeros(n, dtype=np.float64)
    if pairs.size:
        i = pairs[:, 0]
        j = pairs[:, 1]
        counts[i] += 1.0
        counts[j] += 1.0
        di = np.abs(z[i] - z[j])
        # robust per-point scatter of neighbour depth differences (MAD scale)
        acc_abs = np.zeros(n)
        np.add.at(acc_abs, i, di)
        np.add.at(acc_abs, j, di)
        with np.errstate(invalid="ignore", divide="ignore"):
            mean_abs = np.where(counts > 0, acc_abs / np.maximum(counts, 1.0), 0.0)
        scatter = 1.4826 * mean_abs
    tight = counts / (counts + 3.0)
    tol = np.maximum(mad_floor_m, mad_rel * np.abs(z))
    consistent = np.exp(-0.5 * (scatter / tol) ** 2)
    cor = np.where(corroborated, corroborated_boost, 1.0)
    return shallow * persist * tight * consistent * cor


def splat_kde(
    row: np.ndarray,
    col: np.ndarray,
    weight: np.ndarray,
    shape: tuple[int, int],
    *,
    sigma_px: float = 2.0,
    truncate: float = 3.0,
) -> np.ndarray:
    """Bilinear splat of weighted points followed by a truncated Gaussian KDE."""
    grid = np.zeros(shape, dtype=np.float64)
    if row.size:
        r0 = np.floor(row).astype(np.int64)
        c0 = np.floor(col).astype(np.int64)
        ok = (r0 >= 0) & (r0 < shape[0] - 1) & (c0 >= 0) & (c0 < shape[1] - 1)
        r0, c0, wr, wc = r0[ok], c0[ok], row[ok], col[ok]
        w = weight[ok]
        fr = np.clip(wr - r0, 0.0, 1.0)
        fc = np.clip(wc - c0, 0.0, 1.0)
        np.add.at(grid, (r0, c0), w * (1 - fr) * (1 - fc))
        np.add.at(grid, (r0, c0 + 1), w * (1 - fr) * fc)
        np.add.at(grid, (r0 + 1, c0), w * fr * (1 - fc))
        np.add.at(grid, (r0 + 1, c0 + 1), w * fr * fc)
    return ndimage.gaussian_filter(grid, sigma=sigma_px, truncate=truncate, mode="constant", cval=0.0)


def fuse_families(magnetic: np.ndarray, gravity: np.ndarray) -> np.ndarray:
    """Registered fusion: per-family p99.5 normalization then (M+G+sqrt(MG))/3."""
    def norm(field: np.ndarray) -> np.ndarray:
        pos = field[field > 0]
        if pos.size == 0:
            return field
        return np.clip(field / np.quantile(pos, 0.995), 0.0, 1.0)

    m = norm(magnetic)
    g = norm(gravity)
    return (m + g + np.sqrt(m * g)) / 3.0


def catalogue_ramp(labels: np.ndarray, resolution_m: float = RESOLUTION_M,
                   zero_m: float = 150.0, full_m: float = 450.0) -> np.ndarray:
    """Soft flank attenuation: 0 at <= zero_m from known pixels, 1 at >= full_m."""
    known = np.asarray(labels) == 1
    if not known.any():
        return np.ones(labels.shape, dtype=np.float64)
    dist_m = ndimage.distance_transform_edt(~known) * resolution_m
    return np.clip((dist_m - zero_m) / (full_m - zero_m), 0.0, 1.0)


def load_h4_solutions(csv_path, transform_origin: tuple[float, float] = (243350.0, 4508550.0),
                      resolution_m: float = RESOLUTION_M) -> dict[str, np.ndarray]:
    """Read the committed H4 Euler solution cloud; UTM metres -> grid pixels."""
    import gzip
    import csv as _csv

    fam: dict[str, list[np.ndarray]] = {"tmi": [], "gravity": []}
    with gzip.open(csv_path, "rt") as handle:
        reader = _csv.DictReader(handle)
        for rec in reader:
            family = rec["family"]
            key = "tmi" if family == "tmi" else "gravity"
            easting = float(rec["easting_m"])
            northing = float(rec["northing_m"])
            col = (easting - transform_origin[0]) / resolution_m
            row = (transform_origin[1] - northing) / resolution_m
            fam[key].append((row, col, float(rec["effective_depth_m"])))
    out = {}
    for key, rows in fam.items():
        if rows:
            arr = np.asarray(rows, dtype=np.float64)
            out[key] = {"row": arr[:, 0], "col": arr[:, 1], "depth_m": arr[:, 2]}
        else:
            out[key] = {"row": np.empty(0), "col": np.empty(0), "depth_m": np.empty(0)}
    return out


def corroborate(
    row: np.ndarray,
    col: np.ndarray,
    z: np.ndarray,
    solutions: dict[str, np.ndarray],
    *,
    xy_tol_m: float = 300.0,
    depth_floor_m: float = 300.0,
    depth_rel: float = 0.4,
    resolution_m: float = RESOLUTION_M,
) -> np.ndarray:
    """True where an H4 Euler solution of the same family corroborates the peak."""
    flag = np.zeros(row.size, dtype=bool)
    if row.size == 0 or solutions["row"].size == 0:
        return flag
    tree = cKDTree(np.column_stack([solutions["row"], solutions["col"]]) * resolution_m)
    dist, idx = tree.query(np.column_stack([row, col]) * resolution_m, k=1,
                           distance_upper_bound=xy_tol_m)
    found = np.isfinite(dist)
    if found.any():
        tol = np.maximum(depth_floor_m, depth_rel * np.abs(z))
        cand = np.where(found, idx, 0)
        flag = found & (np.abs(solutions["depth_m"][cand] - z) <= tol)
    return flag


def local_depth_stats(z: np.ndarray, sol: np.ndarray, win: int = 5) -> tuple[np.ndarray, np.ndarray]:
    """Per-cell count and std of depth over solution cells in a win x win window."""
    s = sol.astype(np.float64)
    zf = np.where(sol, np.asarray(z, dtype=np.float64), 0.0)
    area = float(win * win)
    cnt = ndimage.uniform_filter(s, size=win, mode="constant", cval=0.0) * area
    m1 = ndimage.uniform_filter(zf, size=win, mode="constant", cval=0.0) * area
    m2 = ndimage.uniform_filter(zf * zf, size=win, mode="constant", cval=0.0) * area
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(cnt > 0, m1 / np.maximum(cnt, 1.0), 0.0)
        var = np.where(cnt > 0, m2 / np.maximum(cnt, 1.0) - mean * mean, 0.0)
    return cnt, np.sqrt(np.maximum(var, 0.0))


def dense_family_weights(
    z: np.ndarray,
    sol: np.ndarray,
    a_cont: np.ndarray,
    z_cont_corrected: np.ndarray,
    corroborated: np.ndarray | None = None,
    *,
    win: int = 5,
    min_window_count: float = 10.0,
    depth_scale_m: float = 1500.0,
    std_floor_m: float = 250.0,
    std_rel: float = 0.35,
    persist_depth_floor_m: float = 400.0,
    persist_depth_rel: float = 0.4,
    non_persistent_factor: float = 0.35,
    corroborated_boost: float = 2.5,
) -> tuple[np.ndarray, np.ndarray]:
    """Registered dense H8 weight: shallow x persistent x density x depth-consistent x corroborated.

    Returns (weight, persist_mask).
    """
    cnt, std = local_depth_stats(z, sol, win=win)
    tol = np.maximum(std_floor_m, std_rel * np.abs(z))
    consist = np.exp(-0.5 * (std / tol) ** 2)
    consist = np.where(cnt >= min_window_count, consist, 0.0)
    zc = np.asarray(z_cont_corrected, dtype=np.float64)
    depth_tol = np.maximum(persist_depth_floor_m, persist_depth_rel * np.abs(z))
    amp_pos = sol & np.isfinite(a_cont) & np.isfinite(zc)
    thresh = float(np.quantile(a_cont[amp_pos], 0.95)) if amp_pos.any() else np.inf
    persist = sol & np.isfinite(zc) & (np.abs(z - zc) <= depth_tol) & (a_cont >= thresh)
    persist_factor = np.where(persist, 1.0, non_persistent_factor)
    density = cnt / (cnt + 3.0)
    shallow = np.exp(-np.clip(z, 0.0, None) / depth_scale_m)
    cor = np.ones(z.shape, dtype=np.float64)
    if corroborated is not None:
        cor = np.where(corroborated, corroborated_boost, 1.0)
    w = shallow * persist_factor * density * consist * cor
    return np.where(sol, w, 0.0), persist
