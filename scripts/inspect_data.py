#!/usr/bin/env python
"""Re-derive every data claim in this repository from the bytes, and hash-pin the result.

Writes ``data/evidence/submission_format.json`` and ``data/evidence/band_inventory.json``.
Nothing in the documentation is allowed to state a number that this script cannot reproduce.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import rasterio

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gems40.layers import BANDS

#: Official package name -> sha256.  The label raster ships as ``existing_faults.tif``
#: in the DrivenData package and is stored as ``labels.tif`` in this mirror; both names
#: are checked so a renamed local copy cannot silently turn into an all-False mask.
PINNED = {
    "training_features.tif": "4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5",
    "existing_faults.tif": "7ba308ccdc4418b31a178f4f1ef21aaa6e152e4028f2f6f64b01f7eb25ae4093",
    "labels.tif": "7ba308ccdc4418b31a178f4f1ef21aaa6e152e4028f2f6f64b01f7eb25ae4093",
    "example_submission.tif": "2176d08e485aa2cd2860ce8df539db4faf4d76163b38a4dd8c30a40454d35cbc",
}


def sha256(path: Path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def main() -> int:
    ddir = Path(sys.argv[1] if len(sys.argv) > 1 else "data")
    out = Path("data/evidence")
    out.mkdir(parents=True, exist_ok=True)
    report = {"files": {}, "pins_match": {}}
    for name, pin in PINNED.items():
        p = ddir / name
        if not p.exists():
            report["files"][name] = dict(present=False)
            report["pins_match"][name] = None
            continue
        with rasterio.open(p) as ds:
            info = dict(present=True, bytes=p.stat().st_size, sha256=sha256(p),
                        shape=(ds.height, ds.width), count=ds.count, dtype=ds.dtypes[0],
                        crs=ds.crs.to_string() if ds.crs else None,
                        transform=[float(v) for v in tuple(ds.transform)[:6]],
                        nodata=str(ds.nodata))
        report["files"][name] = info
        report["pins_match"][name] = (info["sha256"] == pin)

    # band inventory straight from the GeoTIFF band tags
    feats = ddir / "training_features.tif"
    inv = {}
    if feats.exists():
        with rasterio.open(feats) as ds:
            for i in range(1, ds.count + 1):
                t = ds.tags(i)
                inv[i] = dict(band_name=t.get("band_name"), category=t.get("data_category"),
                              description=t.get("description"))
    report["band_inventory"] = inv
    report["band_inventory_matches_layers_py"] = all(
        inv.get(v["band"], {}).get("band_name") == v["band_name"] for v in BANDS.values()
    ) if inv else None

    # footprint and value conventions (measured from the sample submission)
    sample = ddir / "example_submission.tif"
    if sample.exists():
        with rasterio.open(sample) as ds:
            a = ds.read(1)
        foot = np.isfinite(a)
        lab_path = next((ddir / n for n in ("existing_faults.tif", "labels.tif")
                         if (ddir / n).exists()), None)
        labels = None
        if lab_path is not None:
            with rasterio.open(lab_path) as lds:
                raw = lds.read(1)
                nodata = lds.nodata
            labels = raw > 0
            labels_defined = np.isfinite(raw) if nodata is None else (raw != nodata)
        report["footprint"] = dict(
            footprint_pixels=int(foot.sum()), nan_outside=int(np.isnan(a[~foot]).sum()),
            outside_all_nan=bool(np.isnan(a[~foot]).all()),
            inside_unique_values=[float(v) for v in np.unique(a[foot])][:10],
            inside_min=float(a[foot].min()), inside_max=float(a[foot].max()),
            all_values_in_0_1=bool((a[foot] >= 0).all() and (a[foot] <= 1).all()),
            footprint_matches_labels_defined=bool(
                labels is not None and foot.shape == labels_defined.shape
                and np.array_equal(foot, labels_defined)),
            footprint_only_px=int((foot & ~labels_defined).sum()) if labels is not None else None,
            labels_defined_only_px=int((labels_defined & ~foot).sum()) if labels is not None else None,
            labels_source=lab_path.name if lab_path is not None else None,
        )
        if labels is not None:
            report["footprint"]["catalogue_positive_px"] = int(labels.sum())
    (out / "submission_format.json").write_text(json.dumps(report, indent=2))
    (out / "band_inventory.json").write_text(json.dumps(inv, indent=2))
    print(json.dumps(report["pins_match"], indent=2))
    print(json.dumps(report.get("footprint", {}), indent=2))
    print(f"[out] {out/'submission_format.json'}")
    return 0 if all(v is not False for v in report["pins_match"].values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
