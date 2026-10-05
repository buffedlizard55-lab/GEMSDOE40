"""Grid input/output, footprint handling and submission-format validation.

The submission contract (verified against the sample submission's own bytes and the problem page):
  * single-band GeoTIFF, float32
  * shape 3730 x 3292 (rows x cols), EPSG:32611
  * transform (100, 0, 243350, 0, -100, 4508550)  -> 100 m pixels, north-up
  * values in [0, 1] inside the valid footprint
  * NaN outside the valid footprint (the sample submission is exactly 0.0/1.0 inside its
    5,167,373-pixel footprint and NaN in the other 7,111,787 pixels)

Measured, not assumed: ``scripts/inspect_data.py`` re-derives every one of those numbers from
the bytes and writes ``data/evidence/submission_format.json``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio

SHAPE = (3730, 3292)
CRS_EPSG = 32611
TRANSFORM = (100.0, 0.0, 243350.0, 0.0, -100.0, 4508550.0)
NODATA_SENTINEL = -3.4028234663852886e38


def read_band(path: str | Path, band: int = 1, nodata: float | None = NODATA_SENTINEL) -> np.ndarray:
    """Read one band as float64 with NaN for the nodata sentinel."""
    with rasterio.open(path) as ds:
        a = ds.read(band).astype(np.float64)
        nd = ds.nodata if ds.nodata is not None else nodata
    if nd is not None:
        a = np.where(np.isclose(a, nd), np.nan, a)
    return a


def read_profile(path: str | Path) -> dict:
    with rasterio.open(path) as ds:
        return dict(shape=(ds.height, ds.width), count=ds.count, dtype=ds.dtypes[0],
                    crs=ds.crs.to_string() if ds.crs else None, transform=tuple(ds.transform)[:6],
                    nodata=ds.nodata, bounds=tuple(ds.bounds))


def footprint_from_sample(path: str | Path) -> np.ndarray:
    """Valid-footprint mask, taken from the sample submission (finite pixels)."""
    with rasterio.open(path) as ds:
        a = ds.read(1)
    return np.isfinite(a)


def normalise_to_unit(field: np.ndarray, footprint: np.ndarray, floor_quantile: float = 0.0,
                      cap: float = 1.0) -> np.ndarray:
    """Clip negative density to zero, scale the footprint max to ``cap``, NaN outside."""
    f = np.where(footprint, np.nan_to_num(field, nan=0.0, posinf=0.0, neginf=0.0), np.nan)
    f = np.maximum(f, 0.0)
    if floor_quantile > 0.0:
        thr = np.nanquantile(f[footprint], floor_quantile)
        f = np.where(f >= thr, f, 0.0)
    m = np.nanmax(f[footprint])
    if not np.isfinite(m) or m <= 0:
        raise ValueError("field has no positive mass inside the footprint")
    f = f * (cap / m)
    return np.clip(f, 0.0, 1.0)


def write_submission(path: str | Path, values: np.ndarray, footprint: np.ndarray,
                     nan_outside: bool = True, dtype: str = "float32") -> dict:
    """Write the submission raster and return an audit receipt for it."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    v = np.asarray(values, dtype=np.float64).copy()
    if v.shape != SHAPE:
        raise ValueError(f"shape {v.shape} != required {SHAPE}")
    inside = footprint
    if not np.isfinite(v[inside]).all():
        raise ValueError("non-finite value inside the footprint")
    if (v[inside] < 0).any() or (v[inside] > 1).any():
        raise ValueError("value outside [0, 1] inside the footprint")
    outside = ~inside
    v[outside] = np.nan if nan_outside else 0.0
    out = v.astype(dtype)
    profile = dict(driver="GTiff", height=SHAPE[0], width=SHAPE[1], count=1, dtype=dtype,
                   crs=rasterio.crs.CRS.from_epsg(CRS_EPSG), transform=rasterio.Affine(*TRANSFORM),
                   nodata=None, compress="deflate", predictor=3, tiled=False)
    with rasterio.open(path, "w", **profile) as ds:
        ds.write(out, 1)
    return check_submission(path, footprint)


def check_submission(path: str | Path, footprint: np.ndarray) -> dict:
    """The 12 portal checks: geometry, dtype, range, footprint conventions."""
    import hashlib

    p = Path(path)
    sha = hashlib.sha256(p.read_bytes()).hexdigest()
    with rasterio.open(p) as ds:
        a = ds.read(1)
        prof = dict(shape=(ds.height, ds.width), count=ds.count, dtype=ds.dtypes[0],
                    crs=ds.crs.to_string() if ds.crs else None, transform=tuple(ds.transform)[:6],
                    nodata=ds.nodata)
    inside = footprint
    outside = ~footprint
    vals = a[inside]
    checks = {
        "single_band": prof["count"] == 1,
        "dtype_float32": prof["dtype"] == "float32",
        "dimensions_3730x3292": prof["shape"] == SHAPE,
        "crs_epsg_32611": prof["crs"] == "EPSG:32611",
        "transform_exact": tuple(round(x, 6) for x in prof["transform"]) == TRANSFORM,
        "in_footprint_all_finite": bool(np.isfinite(vals).all()),
        "in_footprint_range_0_1": bool((vals >= 0).all() and (vals <= 1).all()),
        "outside_footprint_all_nan": bool(np.isnan(a[outside]).all()),
        "outside_footprint_zero_or_nan": bool((np.isnan(a[outside]) | (a[outside] == 0)).all()),
        "validator_range_0_1_guaranteed": bool(np.isfinite(vals).all() and (vals >= 0).all()
                                               and (vals <= 1).all()),
        "nonzero_support": int((vals > 0).sum()),
        "max_value": float(vals.max()) if vals.size else 0.0,
    }
    receipt = dict(filename=p.name, size_bytes=p.stat().st_size, sha256=sha, **prof,
                   footprint_pixels=int(inside.sum()), emitted_pixels=int((vals > 0).sum()),
                   emitted_mass=float(vals.sum()),
                   in_footprint_min=float(vals.min()) if vals.size else None,
                   in_footprint_max=float(vals.max()) if vals.size else None,
                   nan_outside_footprint=int(np.isnan(a[outside]).sum()),
                   checks=checks,
                   all_checks_passed=bool(all(v for k, v in checks.items()
                                              if k not in ("nonzero_support", "max_value"))))
    return receipt
