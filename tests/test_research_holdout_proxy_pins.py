import hashlib

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

from gemsdoe40.research_holdout import EXPECTED_PROXY_SHA256, read_proxy_truth


def _write(path, array, *, dtype, nodata, transform, crs="EPSG:32611"):
    array = np.asarray(array)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=array.shape[0],
        width=array.shape[1],
        count=1,
        dtype=dtype,
        crs=crs,
        transform=transform,
        nodata=nodata,
    ) as dst:
        dst.write(array.astype(dtype), 1)


def test_current_and_historical_proxy_pins_are_explicitly_separate(tmp_path):
    transform = from_origin(1000, 2000, 100, 100)
    sample = tmp_path / "sample.tif"
    labels = tmp_path / "labels.tif"
    proxy = tmp_path / "proxy.tif"
    _write(sample, [[0.0, 0.0], [np.nan, 0.0]], dtype="float32", nodata=np.nan, transform=transform)
    _write(labels, [[1, 0], [-1, 0]], dtype="int8", nodata=-1, transform=transform)
    _write(proxy, [[1, 1], [0, 0]], dtype="uint8", nodata=None, transform=transform)
    digest = hashlib.sha256(proxy.read_bytes()).hexdigest()
    assert digest != EXPECTED_PROXY_SHA256

    truth, valid, _, info = read_proxy_truth(
        proxy,
        sample,
        labels,
        expected_proxy_sha256=digest,
    )
    assert info["proxy_generation"] == "explicitly_pinned_other"
    assert valid.sum() == 3
    assert truth.sum() == 1  # exact known-label pixel is removed; no distance buffer

    with pytest.raises(ValueError, match="SGMC proxy hash mismatch"):
        read_proxy_truth(proxy, sample, labels)  # default remains the old H2-B pin
