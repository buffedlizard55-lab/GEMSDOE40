"""Official fault-catalogue handling (USGS Quaternary faults + INGENIOUS).

The competition's public label raster is the union of two official sources: the
USGS Quaternary fault and fold database and the INGENIOUS (Initiative for
Geothermal Exploration of New Faults in the US) compilation.  The merged
shapefile published with the competition - ``qfaults_ingenious_nad83conus117`` -
carries per-record provenance text in ``SLIPINFO``, ``RECINFO``, ``GEOCOMM``,
``COMMENTS`` and ``GEOCOMM``, which lets the two provenances be separated
without any private information.

That split is the only *officially sourced* surrogate for the hidden evaluation
set available to a competitor: faults mapped by the INGENIOUS effort are
exactly the kind of structure the organiser holds out, and they sit inside the
same 100 m footprint with the same geometry conventions.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import rasterio
from rasterio import features, warp

# Provenance markers seen in the merged release (measured on the file itself:
# "TRK 2022 ING interp ...", "EJK ING lidar interp ...", "ING Lidar fault mapping
# compilation", ...).  Every record whose provenance text carries one of these
# is treated as INGENIOUS-derived; the rest are classic USGS Quaternary fault
# and fold database records.  ``NAME`` is deliberately not searched: it contains
# ordinary words that can begin or end with "ING".
INGENIOUS_RE = re.compile(r"\bing\b|ingenious", re.IGNORECASE)

TEXT_FIELDS = ("SLIPINFO", "RECINFO", "GEOCOMM", "COMMENTS")


@dataclass
class FaultRecord:
    attrs: dict
    geom: dict          # GeoJSON-like geometry in the source CRS
    source: str = "usgs"          # "usgs" | "ingenious"
    lengths_m: list = field(default_factory=list)


def _text(attrs: dict) -> str:
    return " ".join(str(attrs.get(k) or "") for k in TEXT_FIELDS).lower()


def source_crs(shp_path: str | Path) -> str:
    """WKT CRS of the shapefile (rasterio is built without an OGR driver here)."""
    wkt = Path(shp_path).with_suffix(".prj").read_text()
    return rasterio.crs.CRS.from_wkt(wkt)


def classify_source(attrs: dict) -> str:
    """Return ``"ingenious"`` when the record's provenance text says so."""
    return "ingenious" if INGENIOUS_RE.search(_text(attrs)) else "usgs"


def read_qfaults(shp_path: str | Path) -> list[FaultRecord]:
    """Read the merged QFAULT/INGENIOUS shapefile into records."""
    import shapefile  # pyshp

    reader = shapefile.Reader(str(shp_path))
    fields = [f[0] for f in reader.fields[1:]]
    out: list[FaultRecord] = []
    for sr in reader.iterShapeRecords():
        attrs = dict(zip(fields, sr.record))
        geom = sr.shape.__geo_interface__
        out.append(FaultRecord(attrs=attrs, geom=geom, source=classify_source(attrs)))
    return out


def to_crs(records: list[FaultRecord], src_crs, dst_crs, *, bbox=None) -> list[FaultRecord]:
    """Reproject geometries; ``bbox`` (dst CRS) is applied after reprojection."""
    out = []
    for rec in records:
        try:
            geom = warp.transform_geom(src_crs, dst_crs, rec.geom)
        except Exception:
            continue
        if bbox is not None:
            xs, ys = _coords(geom)
            if xs.size == 0 or xs.max() < bbox[0] or xs.min() > bbox[2] \
                    or ys.max() < bbox[1] or ys.min() > bbox[3]:
                continue
        out.append(FaultRecord(attrs=rec.attrs, geom=geom, source=rec.source))
    return out


def _coords(geom: dict) -> tuple[np.ndarray, np.ndarray]:
    if geom["type"] == "LineString":
        arr = np.asarray(geom["coordinates"], dtype=float)
    elif geom["type"] == "MultiLineString":
        arr = np.concatenate([np.asarray(c, dtype=float) for c in geom["coordinates"]])
    elif geom["type"] == "Polygon":
        arr = np.asarray(geom["coordinates"][0], dtype=float)
    elif geom["type"] == "MultiPolygon":
        arr = np.concatenate([np.asarray(c[0], dtype=float) for c in geom["coordinates"]])
    elif geom["type"] in ("Point", "MultiPoint"):
        arr = np.asarray(geom["coordinates"], dtype=float).reshape(-1, 2)
    else:
        raise ValueError(f"unsupported geometry {geom['type']}")
    return arr[:, 0], arr[:, 1]


def rasterize(records: list[FaultRecord], *, shape, transform, source: str | None = None,
              all_touched: bool = False) -> np.ndarray:
    """Rasterize records onto a grid, optionally only one provenance."""
    sel = [r for r in records if source is None or r.source == source]
    if not sel:
        return np.zeros(shape, dtype=bool)
    return features.rasterize(
        ((r.geom, 1) for r in sel), out_shape=shape, transform=transform, fill=0,
        all_touched=all_touched, dtype="uint8").astype(bool)


def source_report(records: list[FaultRecord]) -> dict:
    """Counts of records and mapped length per provenance."""
    counts = {"usgs": 0, "ingenious": 0}
    lengths = {"usgs": 0.0, "ingenious": 0.0}
    for r in records:
        counts[r.source] += 1
        try:
            xs, ys = _coords(r.geom)
            lengths[r.source] += float(np.hypot(np.diff(xs), np.diff(ys)).sum())
        except ValueError:
            pass
    return dict(records=counts, length_m=lengths)
