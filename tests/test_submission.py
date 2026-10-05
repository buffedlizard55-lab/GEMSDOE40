"""Format contract: [0, 1], float32, EPSG:32611, 3730×3292, no portal-poison values."""
from __future__ import annotations

import numpy as np
from pathlib import Path

import rasterio
from rasterio.transform import Affine

from gemsdoe40.grid import write_submission, check_submission, sha256
from gemsdoe40 import GRID_HEIGHT, GRID_WIDTH, TRANSFORM


def _tiny_template(tmp_path: Path, h=16, w=12):
    path = tmp_path / "template.tif"
    transform = Affine(100.0, 0.0, 243350.0, 0.0, -100.0, 4508550.0)
    data = np.zeros((h, w), dtype=np.float32)
    data[:, : w // 2] = np.nan  # left half outside footprint
    with rasterio.open(
        path, "w", driver="GTiff", height=h, width=w, count=1,
        dtype="float32", crs="EPSG:32611", transform=transform, nodata=np.nan,
    ) as dst:
        dst.write(data, 1)
    footprint = np.isfinite(data)
    return path, footprint


def test_zeros_outside_all_finite_in_range(tmp_path):
    template, footprint = _tiny_template(tmp_path)
    vals = np.zeros(footprint.shape, dtype=np.float32)
    vals[footprint] = 0.37
    out = tmp_path / "sub.tif"
    write_submission(vals, template, out, footprint, outside="zero")
    with rasterio.open(out) as s:
        a = s.read(1)
        assert s.count == 1
        assert a.dtype == np.float32
        assert s.crs.to_epsg() == 32611
        assert s.nodata is None
    assert np.isfinite(a).all()
    assert a.min() >= 0.0 and a.max() <= 1.0
    rec = check_submission(out, footprint)
    assert rec["ok_range"]
    assert rec["nan_inside"] == 0
    assert rec["positive_px"] == int(footprint.sum())


def test_rejects_out_of_range(tmp_path):
    template, footprint = _tiny_template(tmp_path)
    vals = np.zeros(footprint.shape, dtype=np.float32)
    vals[footprint] = 1.2
    try:
        write_submission(vals, template, tmp_path / "bad.tif", footprint, outside="zero")
        raise AssertionError("should have rejected")
    except ValueError as e:
        assert "outside [0, 1]" in str(e)


def test_nan_outside_matches_sample_contract(tmp_path):
    template, footprint = _tiny_template(tmp_path)
    vals = np.zeros(footprint.shape, dtype=np.float32)
    vals[footprint] = 0.5
    out = tmp_path / "nan.tif"
    write_submission(vals, template, out, footprint, outside="nan")
    with rasterio.open(out) as s:
        a = s.read(1)
    assert np.isnan(a[~footprint]).all()
    assert np.isfinite(a[footprint]).all()
