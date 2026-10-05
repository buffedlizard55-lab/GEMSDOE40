"""The committed primary download must be a legal all-finite [0, 1] GeoTIFF."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import rasterio

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = ROOT / "docs" / "downloads" / "GEMSDOE40-submission.tif"


@pytest.mark.skipif(not PRIMARY.exists(), reason="primary TIF not built")
def test_primary_is_portal_safe():
    with rasterio.open(PRIMARY) as s:
        assert s.count == 1
        assert s.dtypes[0] == "float32"
        assert s.crs.to_epsg() == 32611
        assert s.width == 3292 and s.height == 3730
        assert tuple(s.transform)[:6] == (100.0, 0.0, 243350.0, 0.0, -100.0, 4508550.0)
        assert s.nodata is None
        a = s.read(1)
    assert a.dtype == np.float32
    assert np.isfinite(a).all()
    assert float(a.min()) >= 0.0
    assert float(a.max()) <= 1.0
    assert int((a > 0).sum()) > 1000
