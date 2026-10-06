"""H4: offset-aware, rank-adaptive contact Euler solution clouds.

Reid et al. 1990 eq. (2): (r-r0).grad(F) = A, not necessarily zero.
North/east/down coordinates are in metres. Gravity uses dG/dz as a local
single-top-edge approximation; it is NOT an exact finite-step inversion.
All settings are fixed in docs/research/h4-preregistration-20261006.md.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
from scipy import fft, ndimage
from scipy.spatial import cKDTree


@dataclass(frozen=True)
class FieldConfig:
    name: str
    continuation_m: float
    derivative_order: int
    windows: tuple[int, int]
    edge_guard_m: float
    max_depth_m: float


FIELDS = (
    FieldConfig("tmi", 200.0, 0, (15, 25), 1500.0, 5000.0),
    FieldConfig("iso_grav_anom", 500.0, 1, (31, 51), 3000.0, 8000.0),
)
SETTINGS = {
    "si": 0, "fit_contact_offset": True, "resolution_m": 100.0,
    "fft_pad": 128, "stride": 5, "rank2_horizontal_ratio": 0.10,
    "max_design_condition": 1000.0, "max_relative_residual": 0.30,
    "min_effective_depth_m": 100.0, "max_depth_se_fraction": 0.25,
    "depth_se_floor_m": 100.0, "cluster_xy_m": 600.0,
    "cluster_depth_floor_m": 250.0, "cluster_depth_fraction": 0.35,
    "min_other_neighbours": 2, "depth_decay_m": 1800.0,
    "kde_sigma_px": 2.0, "kde_truncate": 4.0, "family_percentile": 99.5,
}
CLOUD_DTYPE = np.dtype([
    ("row", "f8"), ("col", "f8"), ("depth_m", "f8"), ("depth_se_m", "f8"),
    ("residual", "f8"), ("condition", "f8"), ("offset", "f8"),
    ("rank", "u1"), ("window", "u1"), ("weight", "f8"),
    ("neighbours", "i4"), ("depth_mad_m", "f8"),
])


def config_dict() -> dict[str, Any]:
    return {"settings": SETTINGS, "fields": [asdict(c) for c in FIELDS]}


def spectral_gradients(field: np.ndarray, valid: np.ndarray, *, cell_m: float = 100.0,
                       height_m: float = 0.0, derivative_order: int = 0,
                       pad: int = 128) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Consistent harmonic derivatives, with north = negative row direction.

    Interpolation is only a numerical boundary condition. The caller excludes
    real-data gaps and an explicit guard; no interpolated observation is fitted.
    """
    field = np.asarray(field, dtype=np.float32)
    valid = np.asarray(valid, dtype=bool)
    if field.ndim != 2 or valid.shape != field.shape or not valid.any():
        raise ValueError("field and nonempty valid mask must have one 2-D shape")
    if not np.isfinite(field[valid]).all() or cell_m <= 0 or height_m < 0 or pad < 1:
        raise ValueError("invalid field, spacing, continuation height or padding")
    if derivative_order not in (0, 1):
        raise ValueError("H4 accepts scalar magnetics or first-derivative gravity only")
    filled = field.copy()
    filled[valid] -= np.mean(field[valid], dtype=np.float64)
    if not valid.all():
        idx = ndimage.distance_transform_edt(~valid, return_distances=False, return_indices=True)
        filled[~valid] = filled[tuple(idx[:, ~valid])]
        del idx
    padded = np.pad(filled, pad, mode="reflect")
    spectrum = fft.rfft2(padded, workers=2)
    ny, nx = padded.shape
    kx = (2 * np.pi * fft.rfftfreq(nx, cell_m)).astype(np.float32)[None, :]
    kr = (2 * np.pi * fft.fftfreq(ny, cell_m)).astype(np.float32)[:, None]
    wave = np.hypot(kx, kr)
    spectrum *= np.exp(-wave * height_m)
    if derivative_order:
        spectrum *= wave
    crop = (slice(pad, -pad), slice(pad, -pad))
    tx = fft.irfft2(spectrum * (1j * kx), s=padded.shape, workers=2)[crop].copy()
    ty = fft.irfft2(spectrum * (-1j * kr), s=padded.shape, workers=2)[crop].copy()
    tz = fft.irfft2(spectrum * wave, s=padded.shape, workers=2)[crop].copy()
    return tx.astype(np.float32), ty.astype(np.float32), tz.astype(np.float32)


def solve_contact_windows(gradients: np.ndarray, x: np.ndarray, y: np.ndarray,
                          *, horizontal_ratio: float = 0.10) -> dict[str, np.ndarray]:
    """Fit batches of SI=0 contact patches, retaining A and identifiable rank.

    gradients: (batch, observations, 3), components east/north/down.
    x,y: centered observation coordinates, one per observation. Rank-2 windows
    identify cross-strike offset and depth only. The along-strike coordinate is
    explicitly anchored to the window center, never represented as measured.
    """
    d = np.asarray(gradients, dtype=np.float64)
    x, y = np.asarray(x, dtype=np.float64).ravel(), np.asarray(y, dtype=np.float64).ravel()
    if d.ndim != 3 or d.shape[2] != 3 or d.shape[1] != len(x) or x.shape != y.shape or len(x) < 5:
        raise ValueError("expected finite (batch, n>=5, 3) gradients and n coordinates")
    if not (np.isfinite(d).all() and np.isfinite(x).all() and np.isfinite(y).all()):
        raise ValueError("non-finite contact observations")
    if not 0 <= horizontal_ratio < 1:
        raise ValueError("horizontal rank ratio must be in [0,1)")
    b, n, _ = d.shape
    rhs = d[:, :, 0] * x + d[:, :, 1] * y
    rhmean = rhs.mean(axis=1)
    rhs -= rhmean[:, None]
    mean = d.mean(axis=1)
    horizontal_moment = np.einsum("bni,bnj->bij", d[:, :, :2], d[:, :, :2]) / n
    he, hv = np.linalg.eigh(horizontal_moment)
    normal = hv[:, :, -1]
    ratio = np.divide(he[:, 0], he[:, 1], out=np.ones(b), where=he[:, 1] > 0)
    rank = np.where(ratio <= horizontal_ratio, 2, 3).astype(np.uint8)
    centered = d - mean[:, None, :]
    position = np.full((b, 3), np.nan)
    errors = np.full((b, 3), np.nan)
    condition = np.full(b, np.inf)
    residual = np.full(b, np.inf)
    for k in (2, 3):
        group = np.flatnonzero(rank == k)
        if not len(group):
            continue
        if k == 2:
            a = np.stack((np.einsum("bni,bi->bn", centered[group, :, :2], normal[group]),
                          centered[group, :, 2]), axis=-1)
        else:
            a = centered[group]
        scale = np.sqrt(np.mean(a * a, axis=1))
        supported = (scale > np.finfo(np.float64).tiny ** 0.25).all(axis=1)
        scale_safe = np.where(scale > 0, scale, 1.0)
        a /= scale_safe[:, None, :]
        gram = np.einsum("bni,bnj->bij", a, a) / n
        cross = np.einsum("bni,bn->bi", a, rhs[group]) / n
        eigen, vectors = np.linalg.eigh(gram)
        good = supported & (eigen[:, 0] > 1e-12)
        if not good.any():
            continue
        gi = group[good]
        eig, vec, cross = eigen[good], vectors[good], cross[good]
        inv = np.einsum("bij,bj,bkj->bik", vec, 1.0 / eig, vec)
        coeff_scaled = np.einsum("bij,bj->bi", inv, cross)
        coeff = coeff_scaled / scale_safe[good]
        fitted = np.einsum("bni,bi->bn", a[good], coeff_scaled)
        err = rhs[gi] - fitted
        sse = np.sum(err * err, axis=1)
        rhs_ss = np.sum(rhs[gi] * rhs[gi], axis=1)
        relative = np.sqrt(np.divide(sse, rhs_ss, out=np.full(len(gi), np.inf), where=rhs_ss > 0))
        var = sse / (n - k - 1)
        se = np.sqrt(np.maximum(np.diagonal(inv, axis1=1, axis2=2) * var[:, None] / n, 0)) / scale_safe[good]
        if k == 2:
            position[gi, :2] = normal[gi] * coeff[:, :1]
            position[gi, 2] = coeff[:, 1]
            errors[gi, :2] = np.abs(normal[gi]) * se[:, :1]
            errors[gi, 2] = se[:, 1]
        else:
            position[gi] = coeff
            errors[gi] = se
        residual[gi] = relative
        condition[gi] = np.sqrt(eig[:, -1] / eig[:, 0])
    offset = rhmean - np.sum(mean * position, axis=1)
    return {"position": position, "se": errors, "offset": offset, "rank": rank,
            "condition": condition, "residual": residual, "horizontal_eigen_ratio": ratio}


def solution_cloud(tx: np.ndarray, ty: np.ndarray, tz: np.ndarray, valid: np.ndarray,
                   config: FieldConfig, window: int, *, stride: int = 5,
                   cell_m: float = 100.0, batch_size: int = 256) -> tuple[np.ndarray, dict]:
    """Fit all full-resolution patches centered on a stride lattice, in bounded RAM."""
    tx, ty, tz = (np.asarray(a, dtype=np.float32) for a in (tx, ty, tz))
    valid = np.asarray(valid, dtype=bool)
    if tx.ndim != 2 or any(a.shape != tx.shape for a in (ty, tz, valid)):
        raise ValueError("all derivative fields and validity must have identical 2-D shape")
    if window < 3 or window % 2 != 1 or stride < 1 or batch_size < 1 or cell_m <= 0:
        raise ValueError("odd window >=3 and positive stride/batch/spacing required")
    half = window // 2
    safe = ndimage.minimum_filter(valid.astype(np.uint8), size=window, mode="constant", cval=0).astype(bool)
    lattice = np.zeros(valid.shape, dtype=bool)
    lattice[half:valid.shape[0]-half:stride, half:valid.shape[1]-half:stride] = True
    rows, cols = np.nonzero(safe & lattice)
    offsets = np.arange(-half, half + 1)
    dr, dc = np.meshgrid(offsets, offsets, indexing="ij")
    x, y = dc.ravel() * cell_m, -dr.ravel() * cell_m
    views = [np.lib.stride_tricks.sliding_window_view(a, (window, window)) for a in (tx, ty, tz)]
    output = []
    totals = {"windows_tested": int(len(rows)), "accepted": 0, "rank2_tested": 0,
              "rejected_condition_or_rank": 0, "rejected_residual": 0,
              "rejected_depth": 0, "rejected_horizontal": 0, "rejected_depth_se": 0}
    for start in range(0, len(rows), batch_size):
        rr, cc = rows[start:start+batch_size], cols[start:start+batch_size]
        d = np.stack([v[rr-half, cc-half].reshape(len(rr), -1) for v in views], axis=-1)
        fitted = solve_contact_windows(d, x, y)
        pos = fitted["position"]
        depth = pos[:, 2] - config.continuation_m
        cond_ok = np.isfinite(pos).all(axis=1) & (fitted["condition"] <= SETTINGS["max_design_condition"])
        res_ok = fitted["residual"] <= SETTINGS["max_relative_residual"]
        dep_ok = (depth >= SETTINGS["min_effective_depth_m"]) & (depth <= config.max_depth_m)
        xy_ok = np.hypot(pos[:, 0], pos[:, 1]) <= window * cell_m / 2
        se_ok = fitted["se"][:, 2] <= np.maximum(SETTINGS["depth_se_floor_m"], SETTINGS["max_depth_se_fraction"] * depth)
        totals["rank2_tested"] += int(np.count_nonzero(fitted["rank"] == 2))
        for name, check in (("condition_or_rank", cond_ok), ("residual", res_ok), ("depth", dep_ok),
                            ("horizontal", xy_ok), ("depth_se", se_ok)):
            totals[f"rejected_{name}"] += int(np.count_nonzero(~check))
        good = cond_ok & res_ok & dep_ok & xy_ok & se_ok
        sol = np.zeros(int(good.sum()), dtype=CLOUD_DTYPE)
        sol["row"] = rr[good] - pos[good, 1] / cell_m
        sol["col"] = cc[good] + pos[good, 0] / cell_m
        sol["depth_m"] = depth[good]
        sol["depth_se_m"] = fitted["se"][good, 2]
        sol["residual"] = fitted["residual"][good]
        sol["condition"] = fitted["condition"][good]
        sol["rank"] = fitted["rank"][good]
        sol["offset"] = fitted["offset"][good]
        sol["window"] = window
        output.append(sol)
    cloud = np.concatenate(output) if output else np.empty(0, dtype=CLOUD_DTYPE)
    totals.update(accepted=len(cloud), rank2_accepted=int(np.count_nonzero(cloud["rank"] == 2)),
                  window=window, family=config.name,
                  rejection_counts_overlap=True)
    return cloud, totals


def cluster_weights(cloud: np.ndarray, *, cell_m: float = 100.0) -> np.ndarray:
    """Shallow, cross-window, depth-consistent local solution support; no labels."""
    result = cloud.copy()
    if not len(cloud):
        return result
    tree = cKDTree(np.column_stack((cloud["row"], cloud["col"])) * cell_m)
    xy_radius = SETTINGS["cluster_xy_m"]
    for start in range(0, len(cloud), 1024):
        points = tree.data[start:start + 1024]
        neighbours = tree.query_ball_point(points, xy_radius, workers=2)
        for offset, ids in enumerate(neighbours):
            i = start + offset
            ids = np.asarray(ids, dtype=int)
            tolerance = max(SETTINGS["cluster_depth_floor_m"], SETTINGS["cluster_depth_fraction"] * cloud["depth_m"][i])
            ids = ids[(ids != i) & (np.abs(cloud["depth_m"][ids] - cloud["depth_m"][i]) <= tolerance)]
            n = len(ids)
            result["neighbours"][i] = n
            if n < SETTINGS["min_other_neighbours"] or not np.any(cloud["window"][ids] != cloud["window"][i]):
                continue
            median = np.median(cloud["depth_m"][ids])
            mad = 1.4826 * np.median(np.abs(cloud["depth_m"][ids] - median))
            depth, se, res = cloud["depth_m"][i], cloud["depth_se_m"][i], cloud["residual"][i]
            weight = (np.exp(-depth / SETTINGS["depth_decay_m"])
                      * np.exp(-0.5 * (res / SETTINGS["max_relative_residual"]) ** 2)
                      * np.exp(-0.5 * (se / max(depth, 100.0) / 0.25) ** 2)
                      * n / (n + 3) * np.exp(-0.5 * (mad / tolerance) ** 2))
            result["weight"][i] = weight
            result["depth_mad_m"][i] = mad
    return result


def bilinear_splat(rows: np.ndarray, cols: np.ndarray, weights: np.ndarray,
                   shape: tuple[int, int]) -> np.ndarray:
    """Deposit point mass at subpixel cell centers; preserve interior-point mass."""
    rows, cols, weights = (np.asarray(x, dtype=np.float64) for x in (rows, cols, weights))
    if rows.shape != cols.shape or rows.shape != weights.shape or rows.ndim != 1:
        raise ValueError("point coordinates and weights must be matching 1-D arrays")
    if not all(np.isfinite(x).all() for x in (rows, cols, weights)) or (weights < 0).any():
        raise ValueError("point coordinates/weights must be finite; weights nonnegative")
    result = np.zeros(shape, dtype=np.float32)
    r0, c0 = np.floor(rows).astype(int), np.floor(cols).astype(int)
    fr, fc = rows - r0, cols - c0
    for dr, rw in ((0, 1 - fr), (1, fr)):
        for dc, cw in ((0, 1 - fc), (1, fc)):
            rr, cc = r0 + dr, c0 + dc
            good = (rr >= 0) & (rr < shape[0]) & (cc >= 0) & (cc < shape[1])
            np.add.at(result, (rr[good], cc[good]), (weights * rw * cw)[good])
    return result


def family_kde(cloud: np.ndarray, footprint: np.ndarray) -> tuple[np.ndarray, dict]:
    retained = cloud["weight"] > 0
    if not retained.any():
        raise ValueError("no clustered solutions in this family; do not backfill")
    raw = bilinear_splat(cloud["row"][retained], cloud["col"][retained], cloud["weight"][retained], footprint.shape)
    kde = ndimage.gaussian_filter(raw, sigma=SETTINGS["kde_sigma_px"],
                                  truncate=SETTINGS["kde_truncate"], mode="constant", cval=0)
    kde[~footprint] = 0
    support = kde[footprint & (kde > 0)]
    if not len(support):
        raise ValueError("no KDE support in footprint")
    scale = float(np.percentile(support, SETTINGS["family_percentile"]))
    kde = np.clip(kde / scale, 0, 1).astype(np.float32)
    return kde, {"accepted_solutions": len(cloud), "clustered_solutions": int(retained.sum()),
                 "rank2_clustered": int(np.count_nonzero(retained & (cloud["rank"] == 2))),
                 "positive_kde_cells": len(support), "raw_kde_normalizer": scale,
                 "depth_m_p05_p50_p95": np.percentile(cloud["depth_m"][retained], [5, 50, 95]).tolist(),
                 "weight_sum": float(cloud["weight"].sum())}


def combine_families(magnetic: np.ndarray, gravity: np.ndarray, footprint: np.ndarray,
                     known: np.ndarray) -> np.ndarray:
    if magnetic.shape != gravity.shape or magnetic.shape != footprint.shape or magnetic.shape != known.shape:
        raise ValueError("family/mask shape mismatch")
    if not all(np.isfinite(a).all() and (a >= 0).all() and (a <= 1).all() for a in (magnetic, gravity)):
        raise ValueError("normalized KDE fields must be finite in [0,1]")
    field = (magnetic + gravity + np.sqrt(magnetic * gravity)) / 3.0
    field[~footprint | known] = 0.0
    maximum = float(field.max())
    if maximum <= 0:
        raise ValueError("empty combined field")
    return np.clip(field / maximum, 0.0, 1.0).astype(np.float32)
