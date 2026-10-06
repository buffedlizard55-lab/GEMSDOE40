import numpy as np
import rasterio

from gemsdoe40.raster import validate_candidate, write_candidate


def test_writer_matches_sample_nodata_and_grid(tmp_path):
    sample = tmp_path / "sample.tif"
    candidate = tmp_path / "candidate.tif"
    transform = rasterio.transform.from_origin(500000, 4500000, 100, 100)
    arr = np.zeros((20, 30), dtype=np.float32)
    arr[:, :3] = np.nan
    with rasterio.open(sample, "w", driver="GTiff", height=20, width=30, count=1, dtype="float32", crs="EPSG:32611", transform=transform, nodata=np.nan) as ds:
        ds.write(arr, 1)
    pred = np.zeros((20, 30), dtype=np.float32)
    pred[5, 5] = 1.0
    write_candidate(candidate, pred, sample)
    receipt = validate_candidate(candidate, sample)
    assert receipt["valid"]
    assert receipt["candidate"]["in_footprint_nonzero"] == 1
    with rasterio.open(candidate) as ds:
        out = ds.read(1)
        assert np.isnan(out[:, :3]).all()
        assert ds.nodata is not None and np.isnan(ds.nodata)
