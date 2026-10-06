"""Template-aware reading, writing, and portal-format validation."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import rasterio


def band_index_by_name(dataset: rasterio.io.DatasetReader, name: str) -> int:
    matches = [i for i in range(1, dataset.count + 1) if dataset.tags(i).get("band_name") == name]
    if len(matches) != 1:
        available = [dataset.tags(i).get("band_name") for i in range(1, dataset.count + 1)]
        raise ValueError(f"expected exactly one band named {name!r}; found {len(matches)}; available={available}")
    return matches[0]


def canonical_pixel_sha256(array: np.ndarray, outside: np.ndarray) -> str:
    """Hash canonical float32 pixels, normalizing all outside-NaN payloads."""
    values = np.asarray(array, dtype=np.float32).copy()
    outside = np.asarray(outside, dtype=bool)
    if values.shape != outside.shape:
        raise ValueError("array and outside mask shapes differ")
    values[outside] = np.float32(-9999.0)
    values = np.ascontiguousarray(values.astype("<f4", copy=False))
    return hashlib.sha256(values.tobytes()).hexdigest()


def validate_candidate(candidate_path: str | Path, template_path: str | Path) -> dict[str, Any]:
    """Re-open a candidate and verify it against the exact supplied sample grid."""
    candidate_path, template_path = Path(candidate_path), Path(template_path)
    with rasterio.open(template_path) as template:
        template_band = template.read(1)
        footprint = np.isfinite(template_band)
        template_meta = {
            "width": template.width,
            "height": template.height,
            "count": template.count,
            "crs": template.crs.to_string() if template.crs else None,
            "transform": list(template.transform)[:6],
            "dtype": template.dtypes[0],
            "nodata": "NaN" if template.nodata is not None and np.isnan(template.nodata) else template.nodata,
            "valid_pixels": int(footprint.sum()),
            "outside_pixels": int((~footprint).sum()),
        }
        with rasterio.open(candidate_path) as candidate:
            arr = candidate.read(1)
            issues: list[str] = []
            if candidate.count != 1:
                issues.append(f"band_count={candidate.count}, expected 1")
            if arr.shape != template_band.shape:
                issues.append(f"shape={arr.shape}, expected={template_band.shape}")
            if candidate.crs != template.crs:
                issues.append(f"crs={candidate.crs}, expected={template.crs}")
            if candidate.transform != template.transform:
                issues.append(f"transform={candidate.transform}, expected={template.transform}")
            if candidate.dtypes[0] != "float32":
                issues.append(f"dtype={candidate.dtypes[0]}, expected float32")
            if candidate.nodata is None or not np.isnan(candidate.nodata):
                issues.append(f"nodata={candidate.nodata!r}, expected NaN to match sample")
            if arr.shape == footprint.shape:
                if not np.isfinite(arr[footprint]).all():
                    issues.append("non-finite value inside sample footprint")
                if np.isnan(arr[~footprint]).sum() != int((~footprint).sum()):
                    issues.append("outside-footprint cells are not all NaN")
                if np.isinf(arr[~footprint]).any():
                    issues.append("infinite outside-footprint cell")
                if footprint.any():
                    finite_inside = arr[footprint]
                    if finite_inside.size and (float(finite_inside.min()) < 0.0 or float(finite_inside.max()) > 1.0):
                        issues.append(f"in-footprint range=[{float(finite_inside.min())}, {float(finite_inside.max())}], expected [0,1]")
                    if not np.isfinite(finite_inside).all():
                        issues.append("non-finite in-footprint prediction")
            nonzero = int(np.count_nonzero(arr[footprint])) if arr.shape == footprint.shape else None
            sha = hashlib.sha256(candidate_path.read_bytes()).hexdigest()
            pixel_sha = canonical_pixel_sha256(arr, ~footprint) if arr.shape == footprint.shape else None
            return {
                "path": str(candidate_path),
                "bytes": candidate_path.stat().st_size,
                "sha256": sha,
                "canonical_pixels_sha256": pixel_sha,
                "template": template_meta,
                "candidate": {
                    "width": candidate.width,
                    "height": candidate.height,
                    "count": candidate.count,
                    "crs": candidate.crs.to_string() if candidate.crs else None,
                    "transform": list(candidate.transform)[:6],
                    "dtype": candidate.dtypes[0],
                    "nodata": "NaN" if candidate.nodata is not None and np.isnan(candidate.nodata) else candidate.nodata,
                    "in_footprint_min": float(np.min(arr[footprint])) if arr.shape == footprint.shape and footprint.any() else None,
                    "in_footprint_max": float(np.max(arr[footprint])) if arr.shape == footprint.shape and footprint.any() else None,
                    "in_footprint_nonzero": nonzero,
                    "outside_nan": int(np.isnan(arr[~footprint]).sum()) if arr.shape == footprint.shape else None,
                },
                "valid": not issues,
                "issues": issues,
            }


def write_candidate(
    output_path: str | Path,
    prediction: np.ndarray,
    template_path: str | Path,
    *,
    description: str = "GEMSDOE40 preregistered H1 joint magnetic-gravity Euler depth-consensus KDE",
) -> None:
    """Write one float32 band using the sample's grid and NaN nodata semantics."""
    output_path, template_path = Path(output_path), Path(template_path)
    prediction = np.asarray(prediction, dtype=np.float32)
    with rasterio.open(template_path) as template:
        template_band = template.read(1)
        footprint = np.isfinite(template_band)
        if prediction.shape != template_band.shape:
            raise ValueError(f"prediction shape {prediction.shape} does not match template {template_band.shape}")
        if not np.isfinite(prediction[footprint]).all():
            raise ValueError("all predictions inside the template footprint must be finite")
        if prediction[footprint].size and ((prediction[footprint] < 0).any() or (prediction[footprint] > 1).any()):
            raise ValueError("predictions inside the template footprint must be in [0,1]")
        out = prediction.copy()
        out[~footprint] = np.nan
        profile = template.profile.copy()
        profile.update(
            driver="GTiff",
            count=1,
            dtype="float32",
            nodata=np.nan,
            compress="DEFLATE",
            predictor=3,
            BIGTIFF="IF_SAFER",
        )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with rasterio.open(output_path, "w", **profile) as dst:
            dst.write(out.astype(np.float32, copy=False), 1)
            dst.update_tags(
                1,
                description=description,
                value_domain="[0,1]; NaN outside sample footprint",
            )


def write_json(path: str | Path, document: dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
