"""Regression for three legacy GEMSDOE23 predictions marked nodata=0."""
import hashlib
import json
import numpy as np
import rasterio
from gemsdoe40.contact_audit import audit
from tests.test_contact_audit import files


def test_nodata_zero_does_not_hide_raw_absence_predictions(tmp_path):
    sample, candidate, cache = files(tmp_path)
    legacy = tmp_path / "legacy.tif"
    with rasterio.open(candidate) as ds:
        profile = ds.profile.copy(); arr = np.zeros(ds.shape, dtype=np.float32)
    arr[3, 3] = 1
    profile["nodata"] = 0
    with rasterio.open(legacy, "w", **profile) as ds:
        ds.write(arr, 1)
    data = legacy.read_bytes()
    blob = hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()
    (cache / (blob + ".tif")).write_bytes(data)
    inv = tmp_path / "inv.json"
    inv.write_text(json.dumps({"unique_rasters": [{"git_blob_sha": blob,
        "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data), "artifacts": [{"path": "legacy.tif"}]}]}))
    result = audit(candidate, sample, inv, cache)
    row = result["comparisons"][0]
    assert result["completeness_pass"]
    assert row["comparison_coverage"] == 1
    assert row["pixels_compared"] == 132
    assert row["legacy_nodata_zero_compared_as_raw_zero"]
