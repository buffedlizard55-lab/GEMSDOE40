"""Competition grid geometry, footprint, and legal GeoTIFF I/O.

All numbers are measured from the official rasters (byte-verified):

    sample_submission.tif sha256 2176d08e485aa2cd2860ce8df539db4faf4d76163b38a4dd8c30a40454d35cbc
    labels.tif            sha256 7ba308ccdc4418b31a178f4f1ef21aaa6e152e4028f2f6f64b01f7eb25ae4093
    training_features.tif sha256 4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5

    EPSG:32611, 100 m, 3730 × 3292, origin (243350, 4508550)
    footprint 5,167,373 finite cells in the sample submission
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.transform import Affine

from . import GRID_HEIGHT, GRID_WIDTH, NODATA_F32, PIXEL_M, TRANSFORM

SENTINEL = np.float32(NODATA_F32)


def affine() -> Affine:
    a, b, c, d, e, f = TRANSFORM
    return Affine(a, b, c, d, e, f)


def open_template(path: str | Path) -> dict:
    with rasterio.open(path) as s:
        return dict(
            driver="GTiff",
            dtype="float32",
            count=1,
            crs=s.crs,
            transform=s.transform,
            width=s.width,
            height=s.height,
        )


def footprint_from_sample(path: str | Path) -> np.ndarray:
    """Boolean footprint = finite cells of the official sample submission."""
    with rasterio.open(path) as s:
        a = s.read(1)
    return np.isfinite(a)


def footprint_from_features(path: str | Path, band: int = 1) -> np.ndarray:
    with rasterio.open(path) as s:
        a = s.read(band)
    return np.isfinite(a) & (a > SENTINEL * np.float32(0.5))


def read_labels(path: str | Path) -> np.ndarray:
    """Catalogue raster: True on mapped fault (value >= 1)."""
    with rasterio.open(path) as s:
        a = s.read(1)
    return a >= 1


def read_band(path: str | Path, band: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (float32 array with NaN on sentinel, valid mask)."""
    with rasterio.open(path) as s:
        a = s.read(band).astype(np.float32)
    valid = np.isfinite(a) & (a > SENTINEL * np.float32(0.5))
    out = np.where(valid, a, np.nan).astype(np.float32)
    return out, valid


def write_submission(
    values: np.ndarray,
    template_path: str | Path,
    out_path: str | Path,
    footprint: np.ndarray,
    outside: str = "zero",
) -> Path:
    """Write a legal single-band float32 GeoTIFF.

    outside='zero' — all-finite, nodata=None (portal-safe against the
    'Predicted values must be in range [0, 1]' rejection).
    outside='nan'  — NaN outside the footprint, nodata=NaN (sample format).
    """
    meta = open_template(template_path)
    out = np.asarray(values, dtype=np.float32)
    if out.shape != (meta["height"], meta["width"]):
        raise ValueError(f"shape {out.shape} != {(meta['height'], meta['width'])}")
    inside = out[footprint]
    if not np.isfinite(inside).all():
        raise ValueError("non-finite value inside the footprint")
    if (inside < 0).any() or (inside > 1).any():
        raise ValueError("value outside [0, 1] inside the footprint")
    if outside == "nan":
        written = np.where(footprint, out, np.float32(np.nan))
        meta["nodata"] = np.nan
    elif outside == "zero":
        written = np.where(footprint, out, np.float32(0.0))
        meta["nodata"] = None
    else:
        raise ValueError(outside)
    meta["compress"] = "deflate"
    meta["zlevel"] = 9
    meta["tiled"] = True
    meta["predictor"] = 3
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out_path, "w", **meta) as dst:
        dst.write(written, 1)
        dst.crs = CRS.from_epsg(32611)
    return out_path


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_submission(path: str | Path, footprint: np.ndarray) -> dict:
    """Independent re-read of a written file; returns a format receipt."""
    with rasterio.open(path) as s:
        a = s.read(1)
        meta = dict(
            crs=str(s.crs),
            width=s.width,
            height=s.height,
            transform=tuple(s.transform)[:6],
            dtype=str(s.dtypes[0]),
            nodata=s.nodata,
            count=s.count,
        )
    inside = a[footprint]
    outside = a[~footprint]
    finite_in = np.isfinite(inside)
    in_range = finite_in & (inside >= 0) & (inside <= 1)
    return dict(
        path=str(path),
        sha256=sha256(path),
        bytes=Path(path).stat().st_size,
        meta=meta,
        footprint_px=int(footprint.sum()),
        finite_inside=int(finite_in.sum()),
        in_range_inside=int(in_range.sum()),
        nan_inside=int((~np.isfinite(inside)).sum()),
        outside_nan=int(np.isnan(outside).sum()) if np.issubdtype(a.dtype, np.floating) else 0,
        outside_zero=int((outside == 0).sum()),
        positive_px=int(((a > 0) & footprint).sum()),
        min_in=float(np.nanmin(inside)) if finite_in.any() else None,
        max_in=float(np.nanmax(inside)) if finite_in.any() else None,
        ok_shape=bool(a.shape == (GRID_HEIGHT, GRID_WIDTH)),
        ok_range=bool(in_range.sum() == finite_in.sum() == int(footprint.sum())),
        ok_single_band=bool(meta["count"] == 1),
        ok_float32=bool(a.dtype == np.float32),
        pixel_m=PIXEL_M,
    )
