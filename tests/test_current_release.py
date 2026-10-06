import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]


def test_real_download_is_float32_and_bounded_with_exact_declared_footprint():
    """The published H8 GeoTIFF must satisfy the official contract byte for byte."""
    c = json.loads((ROOT / 'docs/data/current-candidate.json').read_text())
    p = ROOT / c['path']
    assert hashlib.sha256(p.read_bytes()).hexdigest() == c['sha256']
    audit = json.loads((ROOT / 'docs/data/h8-lineament-audit.json').read_text())
    assert audit['published']['sha256'] == c['sha256']
    with rasterio.open(p) as ds:
        assert ds.count == 1 and ds.dtypes == ('float32',)
        assert ds.shape == (3730, 3292) and ds.crs.to_epsg() == 32611
        assert tuple(ds.transform)[:6] == (100, 0, 243350, 0, -100, 4508550)
        assert np.isnan(ds.nodata)
        a = ds.read(1)
    finite = np.isfinite(a)
    assert int(finite.sum()) == 5_167_373
    assert int(np.isnan(a).sum()) == 7_111_787
    assert a[finite].min() >= 0.0 and a[finite].max() <= 1.0
    assert np.count_nonzero(a[finite]) == c['dots'] == 40_000
    assert np.unique(a[a > 0]).size == audit['distinct_values'] == 19_923
    # the metric-optimal twin keeps the identical support with every value at 1.0
    twin = c['hard_twin']
    assert twin and (ROOT / twin['path']).is_file()
    assert hashlib.sha256((ROOT / twin['path']).read_bytes()).hexdigest() == twin['sha256']
    with rasterio.open(ROOT / twin['path']) as ds:
        b = ds.read(1)
    assert np.array_equal(np.isfinite(b) & (b > 0), np.isfinite(a) & (a > 0))
    assert set(np.unique(b[np.isfinite(b)])) <= {0.0, 1.0}


def test_real_cloud_record_count_and_depth_labels_match_measured_generation():
    """The downloadable solution cloud is the depth-labelled cloud that built the raster."""
    c = json.loads((ROOT / 'docs/data/current-candidate.json').read_text())
    path = ROOT / c['cloud']['path']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == c['cloud']['sha256']
    receipt = json.loads((ROOT / 'docs/data/h8-lineament-generation.json').read_text())
    expected = (receipt['families']['rtp']['retained']
                + receipt['families']['iso_grav_anom']['retained'])
    count, families = 0, set()
    with gzip.open(path, 'rt') as stream:
        reader = csv.DictReader(stream)
        assert reader.fieldnames[:3] == ['family', 'row', 'col']
        for row in reader:
            count += 1
            families.add(row['family'])
            depth = float(row['depth_m'])
            assert 0.0 < depth <= 5000.0
            assert 0.0 <= float(row['cluster_weight']) <= 1.0
            assert 0.0 <= float(row['lineament_coherence']) <= 1.0
            assert float(row['lineament_weight']) >= 0.0
            assert row['corroborated'] in ('0', '1', 'True', 'False')
    assert count == expected == 124_382
    assert families == {'rtp', 'iso_grav_anom'}


def test_whole_brief_is_embedded_in_readme_without_omissions():
    readme=(ROOT/'README.md').read_text()
    embedded=readme.split('<!-- BEGIN USER BRIEF 20261006 -->')[1].split('<!-- END USER BRIEF 20261006 -->')[0].strip()
    assert embedded == (ROOT/'docs/user-prompt-20261006.md').read_text().strip()
