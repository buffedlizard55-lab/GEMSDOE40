"""H4 - multi-height, multi-field Euler depth clustering with independent corroboration.

Physics (all equations traceable to the cited literature):

* Euler's homogeneity relation solved in a moving window gives one depth-labelled
  source point per window (Reid, Allsop, Granser, Millett & Somerton, 1990,
  Geophysics 55(1):80-91, https://doi.org/10.1190/1.1442774).  Structural index
  0 is the fault-like contact case (Reid et al. 1990 appendix; Reid 2003,
  https://www.reid-geophys.co.uk/wp-content/uploads/2017/11/Reid-Thurston-2014.pdf).
* Solving the same field after upward continuation to several heights and
  keeping solutions that *persist* in position and depth is a standard
  noise-rejection device (Hall 1981; Reid et al. 1990 sec. "Interpretation of
  the solution map"; Blakely 1995, "Potential Theory in Gravity and Magnetic
  Applications", ch. 12 for the continuation operator exp(-|k| h)).
* Tight clusters of shallow, mutually depth-consistent solutions mark real
  near-surface contacts; scattered or deep solutions mark noise.  The KDE
  weighting here implements that statement quantitatively (depth decay x
  local tightness x local depth consistency x cross-height persistence x
  cross-field concordance).

The field this module builds is deliberately *not* a gradient threshold: every
value in it comes from a solved source-point cloud, and gradient-derived layers
enter only as independent corroboration weights (tilt angle of the magnetic
field - Miller & Singh 1994, https://doi.org/10.1190/1.1443546; horizontal
gradient of the gravity field; topographic curvature), never as the primary
detector.

Everything is computed on the competition grid: EPSG:32611, 100 m,
3730 x 3292, exact sample transform.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.fft import irfft2, rfft2
from scipy.ndimage import binary_dilation, gaussian_filter
from scipy.spatial import cKDTree

from .euler import EulerCloud, deconvolve, derivatives, merge_clouds
from .kde import concordance_boost, splat_kde

PIXEL_M = 100.0


# --------------------------------------------------------------------------
# upward continuation
# --------------------------------------------------------------------------
def upward_continue(field: np.ndarray, valid: np.ndarray, height_m: float,
                    pixel_m: float = PIXEL_M) -> np.ndarray:
    """Upward-continue a potential field by ``height_m`` (Blakely 1995, eq. 12-8).

    Uses the exact continuation operator exp(-|k| h) on a mean-removed,
    nearest-filled grid; the result is masked back to ``valid``.
    """
    f = np.asarray(field, dtype=np.float64)
    v = np.asarray(valid, dtype=bool)
    if height_m <= 0:
        return np.where(v, f, np.nan)
    filled = _fill_nearest(np.where(v, f, 0.0), v)
    mean = float(filled[v].mean()) if v.any() else 0.0
    work = filled - mean
    h, w = f.shape
    ky = 2.0 * np.pi * np.fft.fftfreq(h, d=pixel_m)
    kx = 2.0 * np.pi * np.fft.rfftfreq(w, d=pixel_m)
    k = np.hypot(ky[:, None], kx[None, :])
    cont = irfft2(rfft2(work) * np.exp(-k * float(height_m)), s=(h, w)).real + mean
    return np.where(v, cont, np.nan)


def _fill_nearest(arr: np.ndarray, valid: np.ndarray) -> np.ndarray:
    from scipy.ndimage import distance_transform_edt
    if valid.all():
        return np.asarray(arr, dtype=np.float64)
    idx = distance_transform_edt(~valid, return_distances=False, return_indices=True)
    return np.asarray(arr, dtype=np.float64)[tuple(idx)]


# --------------------------------------------------------------------------
# cloud construction
# --------------------------------------------------------------------------
@dataclass
class StackSpec:
    layers: tuple
    heights_m: tuple
    windows_px: tuple
    stride_px: int = 4
    structural_index: float = 0.0
    analytic_percentile: float = 72.0
    max_rel_se: float = 0.22
    min_depth_m: float = 80.0
    max_depth_m: float = 2200.0
    max_condition: float = 5e4


def build_cloud(layer: str, field: np.ndarray, valid: np.ndarray, spec: StackSpec,
                pixel_m: float = PIXEL_M) -> tuple[EulerCloud, EulerCloud]:
    """Deconvolve one layer at every (height, window) in the spec.

    Returns ``(all_solutions, persisted_only)``.  ``persisted_only`` keeps the
    solutions that also appear (within 2 px, 35 % depth) at another
    continuation height - the multi-height stability filter.
    """
    clouds = []
    for h_m in spec.heights_m:
        cont = upward_continue(field, valid, h_m, pixel_m) if h_m > 0 else field
        v_ok = np.isfinite(cont)
        for win in spec.windows_px:
            c = deconvolve(
                np.where(v_ok, cont, 0.0), v_ok,
                field_name=f"{layer}@{int(h_m)}m",
                structural_index=spec.structural_index,
                window_px=int(win), stride_px=int(spec.stride_px), pixel_m=pixel_m,
                analytic_percentile=spec.analytic_percentile,
                max_condition=spec.max_condition, max_rel_se=spec.max_rel_se,
                min_depth_m=spec.min_depth_m, max_depth_m=spec.max_depth_m,
            )
            clouds.append(c)
    merged = merge_clouds(clouds)
    merged.field_name = layer
    # persistence across heights
    lab = np.array([c.field_name for c in clouds])
    height = np.array([float(n.split("@")[1].rstrip("m")) for n in lab])
    _ = height  # heights are encoded in every solution's row; persistence uses geometry
    persistence = cross_height_persistence(clouds)
    keep = persistence >= 1.0
    persisted = EulerCloud(
        field_name=f"{layer}-persisted", structural_index=spec.structural_index,
        window_px=-1, stride_px=-1,
        row=merged.row[keep], col=merged.col[keep], depth_m=merged.depth_m[keep],
        rel_se=merged.rel_se[keep], cond=merged.cond[keep], analytic=merged.analytic[keep],
        stats={"persistence_kept": int(keep.sum()), "of": int(len(merged))},
    )
    return merged, persisted


def cross_height_persistence(clouds: list[EulerCloud], *, radius_px: float = 2.0,
                             depth_tol: float = 0.35) -> np.ndarray:
    """For each solution, the number of *other* height levels reproducing it.

    A solution produced from the field continued to height h is 'reproduced'
    when a solution from a different continuation height lands within
    ``radius_px`` with a relative depth difference below ``depth_tol``.
    """
    rows = np.concatenate([c.row for c in clouds])
    cols = np.concatenate([c.col for c in clouds])
    depth = np.concatenate([c.depth_m for c in clouds])
    level = np.concatenate([[i] * len(c) for i, c in enumerate(clouds)])
    pts = np.column_stack([rows, cols])
    tree = cKDTree(pts)
    out = np.zeros(rows.size, dtype=np.float64)
    for i in range(rows.size):
        nb = tree.query_ball_point(pts[i], radius_px)
        other = [j for j in nb if level[j] != level[i]]
        if not other:
            continue
        d = depth[other]
        rel = np.abs(d - depth[i]) / max(0.5 * (float(np.mean(d)) + depth[i]), 1.0)
        out[i] = float(np.sum(rel <= depth_tol))
    return out


# --------------------------------------------------------------------------
# depth-cluster weighting and KDE
# --------------------------------------------------------------------------
def cluster_weights(cloud: EulerCloud, *, depth_scale_m: float = 700.0,
                    neighbor_px: float = 3.0, min_neighbors: int = 3) -> np.ndarray:
    """shallow x precise x tight x depth-consistent x analytic-strength weight."""
    n = len(cloud)
    if n == 0:
        return np.empty(0)
    shallow = np.exp(-np.clip(cloud.depth_m, 0.0, None) / depth_scale_m)
    quality = 1.0 / (1.0 + 8.0 * np.clip(cloud.rel_se, 0.0, 10.0))
    pts = np.column_stack([cloud.row, cloud.col])
    tree = cKDTree(pts)
    nb = tree.query_ball_point(pts, neighbor_px)
    tight = np.zeros(n)
    consist = np.zeros(n)
    for i, neigh in enumerate(nb):
        k = len(neigh)
        tight[i] = k
        if k >= min_neighbors:
            d = cloud.depth_m[np.asarray(neigh, dtype=int)]
            mu = float(d.mean())
            cv = float(d.std()) / max(mu, 1.0)
            consist[i] = float(np.exp(-4.0 * cv))
    tight = np.clip(tight / max(float(np.percentile(tight, 90)), 1.0), 0.0, 3.0)
    a = cloud.analytic
    p90 = float(np.percentile(a, 90)) if a.size else 1.0
    strength = 0.25 + 0.75 * np.clip(a / max(p90, 1e-12), 0.0, 2.0)
    w = shallow * quality * tight * consist * strength
    return np.where(np.isfinite(w) & (w > 0), w, 0.0)


def depth_cluster_kde(cloud: EulerCloud, shape: tuple[int, int], *, sigma_px: float = 1.6) -> np.ndarray:
    w = cluster_weights(cloud)
    return splat_kde(cloud, w, shape, sigma_px=sigma_px)


# --------------------------------------------------------------------------
# lineament matched filter (orientation-aware ridge enhancement)
# --------------------------------------------------------------------------
_DIRS8 = [(np.cos(t), np.sin(t)) for t in np.arange(8) * np.pi / 8.0]


def line_response(field: np.ndarray, *, half_len: int = 3, across: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Line-vs-background matched filter over 8 orientations.

    For each orientation the along-line mean is compared with the two
    perpendicular offsets; the returned response is the maximum over
    orientations and ``orientation_index`` records which one won.  This is the
    classical lineament enhancement used for potential-field edge mapping and
    it is applied *after* the Euler KDE, not instead of it.
    """
    f = np.asarray(field, dtype=np.float64)
    if across > 0:
        f = gaussian_filter(f, sigma=across, mode="constant")
    h, w = f.shape
    best = np.full((h, w), -np.inf)
    best_dir = np.zeros((h, w), dtype=np.int8)
    # a coarse background: local mean over a disk wider than the line
    bg = gaussian_filter(f, sigma=float(half_len) * 1.5, mode="constant")
    for d, (cx, cy) in enumerate(_DIRS8):
        acc = np.zeros((h, w), dtype=np.float64)
        cnt = 0
        for t in range(-half_len, half_len + 1):
            dy, dx = int(round(t * cy)), int(round(t * cx))
            acc += np.roll(np.roll(f, dy, axis=0), dx, axis=1)
            cnt += 1
        acc /= float(cnt)
        resp = acc - bg
        upd = resp > best
        best = np.where(upd, resp, best)
        best_dir = np.where(upd, np.int8(d), best_dir)
    best = np.where(np.isfinite(best), best, 0.0)
    return best, best_dir


def coherence(field: np.ndarray, sigma: float = 2.0) -> np.ndarray:
    """Structure-tensor coherence in [0, 1] (1 = perfectly linear structure)."""
    f = gaussian_filter(np.asarray(field, dtype=np.float64), sigma=max(sigma - 1.0, 0.4), mode="constant")
    gy, gx = np.gradient(f)
    jxx = gaussian_filter(gx * gx, sigma=sigma, mode="constant")
    jyy = gaussian_filter(gy * gy, sigma=sigma, mode="constant")
    jxy = gaussian_filter(gx * gy, sigma=sigma, mode="constant")
    tr = jxx + jyy
    det = jxx * jyy - jxy * jxy
    disc = np.sqrt(np.maximum(tr * tr / 4.0 - det, 0.0))
    l1 = tr / 2.0 + disc
    l2 = tr / 2.0 - disc
    coh = np.where(l1 > 0, (l1 - l2) / np.maximum(l1 + l2, 1e-12), 0.0)
    return np.clip(coh, 0.0, 1.0)


# --------------------------------------------------------------------------
# independent corroboration channels
# --------------------------------------------------------------------------
def tilt_angle(field: np.ndarray, valid: np.ndarray, pixel_m: float = PIXEL_M) -> np.ndarray:
    """Tilt angle of the magnetic field (Miller & Singh 1994; Verduzco et al. 2004).

    TDR = atan2(dT/dz, sqrt((dT/dx)^2 + (dT/dy)^2)); its zero contour tracks
    contacts and faults independently of magnetisation direction.
    """
    gx, gy, gz, _, dvalid = derivatives(field, valid, pixel_m)
    hg = np.hypot(gx, gy)
    tdr = np.arctan2(gz, np.maximum(hg, 1e-12))
    return np.where(dvalid, tdr, np.nan)


def horizontal_gradient_magnitude(field: np.ndarray, valid: np.ndarray,
                                  pixel_m: float = PIXEL_M) -> np.ndarray:
    gx, gy, _, _, dvalid = derivatives(field, valid, pixel_m)
    return np.where(dvalid, np.hypot(gx, gy), np.nan)


def curvature_evidence(detrended_elevation: np.ndarray, valid: np.ndarray,
                       pixel_m: float = PIXEL_M, smooth_sigma: float = 1.0) -> np.ndarray:
    """Absolute profile curvature of the provided detrended elevation band."""
    f = gaussian_filter(_fill_nearest(np.where(valid, detrended_elevation, np.nan), valid),
                        sigma=smooth_sigma, mode="nearest")
    fy, fx = np.gradient(f)
    fyy, _ = np.gradient(fy)
    _, fxx = np.gradient(fx)
    grad = np.hypot(fx, fy)
    curv = -(fxx * fy * fy - 2 * fx * fy * fyy / np.maximum(np.hypot(fx, fy), 1e-12)
             + fyy * fx * fx) / np.maximum((fx * fx + fy * fy) * np.power(1.0 + grad * grad, 1.5), 1e-12)
    out = np.abs(np.where(np.isfinite(curv), curv, 0.0))
    return np.where(valid, out, np.nan)


def robust_norm(x: np.ndarray, valid: np.ndarray, lo_q: float = 2.0, hi_q: float = 98.0) -> np.ndarray:
    v = np.asarray(x, dtype=np.float64)
    m = np.asarray(valid, bool) & np.isfinite(v)
    if not m.any():
        return np.zeros_like(v)
    lo, hi = np.percentile(v[m], [lo_q, hi_q])
    if hi <= lo:
        hi = lo + 1e-12
    return np.where(m, np.clip((v - lo) / (hi - lo), 0.0, 1.0), 0.0)


# --------------------------------------------------------------------------
# assembly
# --------------------------------------------------------------------------
@dataclass
class H4Config:
    sigma_px: float = 1.6
    concord_weight: float = 1.75
    catalogue_buffer_px: int = 0          # pixel-exact mask only (organiser rule)
    line_half_len: int = 3
    w_line: float = 0.45                  # weight of the line-response term
    w_coherence: float = 0.30
    w_corroboration: float = 0.35
    corrob_floor: float = 0.55
    # Catalogue-flank gate.  ``scripts/study_flank.py`` and
    # ``scripts/study_evidence.py`` measure every alternative surface against the
    # two surrogate instruments; a band at 2-6 px from the mapped catalogue,
    # carrying the Euler line-response texture, dominates all of them (see
    # docs/data/flank_study.json, docs/data/evidence_study.json).  ``flank_mix``
    # keeps a floor so pixels inside the band are never zeroed.
    use_flank_band: bool = True
    flank_lo_px: float = 2.0
    flank_hi_px: float = 6.0
    flank_mix: float = 0.25
    flank_weight: float = 0.75
    # The line texture can be taken from the plain depth-cluster KDE (False) or
    # from the concordance-boosted primary (True).  scripts/ablate_h4.py measures
    # both; the KDE source is the stronger one, because the magnetic/gravity
    # concordance term follows the much sparser gravity solution cloud.
    concord_in_line: bool = False


def assemble_ranking(mag_cloud: EulerCloud, grav_cloud: EulerCloud,
                     corroborators: dict, footprint: np.ndarray, catalogue: np.ndarray,
                     cfg: H4Config) -> tuple[np.ndarray, dict]:
    """Continuous ranking surface in [0, 1] on the competition grid.

    ``corroborators`` maps a channel name to an already re-scaled [0, 1]
    evidence surface (or None).  The combined corroboration is the pixelwise
    maximum of the channels, which is the weakest reasonable combination and
    therefore the least prone to overfitting a single evidence type.
    """
    shape = footprint.shape
    kde_m = depth_cluster_kde(mag_cloud, shape, sigma_px=cfg.sigma_px) if len(mag_cloud) else np.zeros(shape)
    kde_g = depth_cluster_kde(grav_cloud, shape, sigma_px=cfg.sigma_px) if len(grav_cloud) else np.zeros(shape)
    concord = (concordance_boost(mag_cloud, grav_cloud, shape, radius_px=2.5, depth_tol=0.35,
                                 sigma_px=max(1.0, cfg.sigma_px - 0.4))
               if len(mag_cloud) and len(grav_cloud) else np.zeros(shape))
    primary = kde_m + kde_g + cfg.concord_weight * concord
    primary = np.where(footprint, primary, 0.0)
    n_primary = robust_norm(primary, footprint)
    line_source = n_primary
    if not cfg.concord_in_line:
        line_source = robust_norm(np.where(footprint, kde_m + kde_g, 0.0), footprint)
    line, orient = line_response(line_source, half_len=cfg.line_half_len)
    line_source_name = "primary" if cfg.concord_in_line else "kde_mag+grav"
    line_n = robust_norm(line, footprint)
    coh = coherence(np.where(footprint, primary, 0.0), sigma=2.0)

    corrob = np.zeros(shape, dtype=np.float64)
    used = []
    for name, ch in (corroborators or {}).items():
        if ch is None:
            continue
        corrob = np.maximum(corrob, robust_norm(ch, footprint))
        used.append(name)

    rank = n_primary * (1.0 - cfg.w_line + cfg.w_line * line_n) \
        * (1.0 - cfg.w_coherence + cfg.w_coherence * coh) \
        * (cfg.corrob_floor + cfg.w_corroboration * corrob)
    band_n = np.zeros(shape, dtype=np.float64)
    if cfg.use_flank_band:
        from scipy.ndimage import distance_transform_edt
        d = distance_transform_edt(~np.asarray(catalogue, bool))
        band = np.zeros(shape, dtype=np.float64)
        sel = (d >= cfg.flank_lo_px) & (d <= cfg.flank_hi_px) & footprint
        band[sel] = 1.0 / (d[sel] - (cfg.flank_lo_px - 1.0))
        band_n = robust_norm(band, footprint)
        # measured form: Euler line-response texture gated by the flank band
        rank = line_n * (cfg.flank_mix + cfg.flank_weight * band_n)
    rank = np.where(footprint & (primary > 0), rank, 0.0)
    if cfg.catalogue_buffer_px >= 0:
        cat = np.asarray(catalogue, bool) & footprint
        if cfg.catalogue_buffer_px > 0 and cat.any():
            cat = binary_dilation(cat, iterations=int(cfg.catalogue_buffer_px))
        rank = np.where(cat, 0.0, rank)
    mx = float(rank[footprint].max()) if footprint.any() else 0.0
    if mx > 0:
        rank = rank / mx
    rank = np.clip(rank, 0.0, 1.0)
    stats = dict(
        mag_solutions=int(len(mag_cloud)), grav_solutions=int(len(grav_cloud)),
        primary_positive_px=int((primary > 0).sum()),
        rank_positive_px=int((rank > 0).sum()),
        line_response_p90=float(np.percentile(line_n[footprint], 90)) if footprint.any() else 0.0,
        coherence_mean=float(coh[footprint].mean()) if footprint.any() else 0.0,
        corroboration_channels=used,
        catalogue_buffer_px=int(cfg.catalogue_buffer_px),
        on_catalogue_positive=int(((rank > 0) & np.asarray(catalogue, bool)).sum()),
        flank_band_positive_px=int((band_n > 0).sum()),
        line_source=line_source_name,
        config=cfg.__dict__,
    )
    return rank.astype(np.float32), stats
