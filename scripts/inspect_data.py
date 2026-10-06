#!/usr/bin/env python
"""Re-derive dataset and submission-format facts from the locally hash-pinned bytes.

Writes ``data/evidence/submission_format.json`` and ``data/evidence/band_inventory.json``. This
script accepts the actual local canonical names (labels.tif and sample_submission.tif) as well as
the names used in historical manifests; a missing legacy alias no longer silently masks a missing
canonical input.
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
from gems40.pins import PINNED_FILES


INPUTS = {
    "training_features.tif": {
        "aliases": ["training_features.tif"],
        "pin": PINNED_FILES["training_features.tif"]["sha256"],
        "role": PINNED_FILES["training_features.tif"]["note"],
    },
    "labels.tif": {
        "aliases": ["labels.tif", "existing_faults.tif"],
        "pin": PINNED_FILES["existing_faults.tif"]["sha256"],
        "role": PINNED_FILES["existing_faults.tif"]["note"],
    },
    "sample_submission.tif": {
        "aliases": ["sample_submission.tif", "example_submission.tif"],
        "pin": PINNED_FILES["example_submission.tif"]["sha256"],
        "role": PINNED_FILES["example_submission.tif"]["note"],
    },
}
PROXIES = {
    "current_643cbe": {
        "path": "external/derived_sgmc_faults_100m_u8.tif",
        "sha256": PINNED_FILES["derived_sgmc_faults_100m_u8.tif"]["sha256"],
        "role": "current owner-derived SGMC proxy only; not organizer truth",
    },
    "historical_h2b_26d142": {
        "path": "external/archive/derived_sgmc_faults_100m_u8_h2b_26d142.tif",
        "sha256": "26d142c4c93282cd94f6950ab96f22aeff59fbbea523d43d662e76fa1b161b5c",
        "role": "historical H2-B SGMC proxy instrument; distinct byte version; not organizer truth",
    },
}


def sha256(path: Path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def raster_info(path: Path) -> dict:
    with rasterio.open(path) as ds:
        return dict(
            path=str(path),
            present=True,
            bytes=path.stat().st_size,
            sha256=sha256(path),
            shape=[ds.height, ds.width],
            count=ds.count,
            dtype=ds.dtypes[0],
            crs=ds.crs.to_string() if ds.crs else None,
            transform=[float(v) for v in tuple(ds.transform)[:6]],
            nodata="NaN" if ds.nodata is not None and np.isnan(ds.nodata) else ds.nodata,
        )


def main() -> int:
    data_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "data")
    output_dir = Path("data/evidence")
    output_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "inspection_date_utc": "2026-10-06",
        "acquisition_note": "Inputs are byte-verified public sibling-repository mirrors; the official DrivenData data tab requires login in this environment.",
        "files": {},
        "pins_match": {},
    }
    resolved: dict[str, Path] = {}
    for canonical, spec in INPUTS.items():
        path = next((data_dir / name for name in spec["aliases"] if (data_dir / name).exists()), None)
        if path is None:
            report["files"][canonical] = {"present": False, "aliases_checked": spec["aliases"]}
            report["pins_match"][canonical] = False
            continue
        resolved[canonical] = path
        info = raster_info(path)
        info["role"] = spec["role"]
        report["files"][canonical] = info
        report["pins_match"][canonical] = info["sha256"] == spec["pin"]

    for label, spec in PROXIES.items():
        path = data_dir / spec["path"]
        if not path.exists():
            report.setdefault("proxies", {})[label] = {"present": False, "path": str(path), "expected_sha256": spec["sha256"], "role": spec["role"]}
            continue
        info = raster_info(path)
        info["expected_sha256"] = spec["sha256"]
        info["pin_matches"] = info["sha256"] == spec["sha256"]
        info["role"] = spec["role"]
        report.setdefault("proxies", {})[label] = info

    band_inventory: dict[int, dict] = {}
    feature_path = resolved.get("training_features.tif")
    if feature_path is not None:
        with rasterio.open(feature_path) as ds:
            for i in range(1, ds.count + 1):
                tags = ds.tags(i)
                band_inventory[i] = {
                    "band_name": tags.get("band_name"),
                    "category": tags.get("data_category"),
                    "description": tags.get("description"),
                }
    report["band_inventory"] = band_inventory
    report["band_inventory_matches_layers_py"] = all(
        band_inventory.get(spec["band"], {}).get("band_name") == spec["band_name"]
        for spec in BANDS.values()
    ) if band_inventory else None

    sample_path = resolved.get("sample_submission.tif")
    labels_path = resolved.get("labels.tif")
    if sample_path is not None:
        with rasterio.open(sample_path) as sample_ds:
            sample = sample_ds.read(1)
            footprint = np.isfinite(sample)
            report["sample_grid"] = {
                "height": sample_ds.height,
                "width": sample_ds.width,
                "crs": sample_ds.crs.to_string() if sample_ds.crs else None,
                "transform": [float(v) for v in tuple(sample_ds.transform)[:6]],
                "resolution_m": float(abs(sample_ds.transform.a)),
                "footprint_pixels": int(footprint.sum()),
                "outside_pixels": int((~footprint).sum()),
                "outside_nan": int(np.isnan(sample[~footprint]).sum()),
                "inside_min": float(sample[footprint].min()),
                "inside_max": float(sample[footprint].max()),
            }
        if labels_path is not None:
            with rasterio.open(labels_path) as ds:
                labels = ds.read(1)
            report["sample_grid"]["known_label_positive_pixels"] = int(np.count_nonzero((labels == 1) & footprint))
            report["sample_grid"]["known_label_exact_overlap_outside"] = int(np.count_nonzero((labels == 1) & ~footprint))
            report["sample_grid"]["labels_aligned"] = bool(labels.shape == footprint.shape)

    (output_dir / "submission_format.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output_dir / "band_inventory.json").write_text(json.dumps(band_inventory, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"pins_match": report["pins_match"], "proxy_pins": {k: v.get("pin_matches") for k, v in report.get("proxies", {}).items()}, "sample_grid": report.get("sample_grid")}, indent=2))
    required = all(report["pins_match"].values()) and bool(report.get("band_inventory_matches_layers_py"))
    required = required and all(item.get("pin_matches") is True for item in report.get("proxies", {}).values())
    print(f"[out] {output_dir / 'submission_format.json'}")
    return 0 if required else 1


if __name__ == "__main__":
    raise SystemExit(main())
