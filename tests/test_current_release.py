"""Current-candidate (H8) release assertions: exact disk re-read of the published artifact."""
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]


def test_real_download_is_float32_and_bounded_with_exact_declared_footprint():
    c = json.loads((ROOT / 'docs/data/current-candidate.json').read_text())
    p = ROOT / c['path']
    assert hashlib.sha256(p.read_bytes()).hexdigest() == c['sha256']
    with rasterio.open(p) as ds:
        assert ds.count == 1 and ds.dtypes == ('float32',)
        assert ds.shape == (3730, 3292) and ds.crs.to_epsg() == 32611
        assert tuple(ds.transform)[:6] == (100, 0, 243350, 0, -100, 4508550)
        assert np.isnan(ds.nodata)
        a = ds.read(1)
    finite = np.isfinite(a)
    assert int(finite.sum()) == 5_167_373
    assert int(np.isnan(a).sum()) == 7_111_787
    assert a[finite].min() == 0 and a[finite].max() == 1
    # 741,536 positive cells; 724,844 distinct float32 values (the float64 field
    # had 731,858 distinct positive values before float32 quantization) plus 0.0
    assert np.count_nonzero(a[finite]) == 741_536
    assert np.unique(a[finite]).size == 724_844


def test_portal_safe_twin_is_all_finite_without_nodata_tag():
    c = json.loads((ROOT / 'docs/data/current-candidate.json').read_text())
    p = ROOT / 'docs' / 'downloads' / c['portal_safe_twin']['filename']
    assert hashlib.sha256(p.read_bytes()).hexdigest() == c['portal_safe_twin']['sha256']
    with rasterio.open(p) as ds:
        assert ds.nodata is None
        a = ds.read(1)
    assert np.isfinite(a).all()
    assert float(a.min()) >= 0.0 and float(a.max()) <= 1.0
    assert np.count_nonzero(a) == 741_536


def test_zip_contains_exactly_one_geotiff_matching_the_zeros_twin():
    import zipfile

    c = json.loads((ROOT / 'docs/data/current-candidate.json').read_text())
    zp = ROOT / 'docs' / 'downloads' / c['zip']['filename']
    assert hashlib.sha256(zp.read_bytes()).hexdigest() == c['zip']['sha256']
    with zipfile.ZipFile(zp) as archive:
        names = archive.namelist()
        assert len(names) == 1 and names[0].endswith('.tif')
        inner = archive.read(names[0])
    zeros = (ROOT / 'docs' / 'downloads' / c['portal_safe_twin']['filename']).read_bytes()
    assert hashlib.sha256(inner).hexdigest() == hashlib.sha256(zeros).hexdigest()


def test_real_cloud_record_count_and_depth_labels_match_measured_generation():
    c = json.loads((ROOT / 'docs/data/current-candidate.json').read_text())
    path = ROOT / c['cloud']['path']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == c['cloud']['sha256']
    count = 0
    fams = {'magnetic': 0, 'gravity': 0}
    with gzip.open(path, 'rt') as stream:
        for row in csv.DictReader(stream):
            count += 1
            family = row['family']
            assert family in fams
            depth = float(row['effective_depth_m'])
            assert 100 <= depth <= 8000
            # raw pre-KDE cell weights: shallowness*persistence*density*consistency*corroboration,
            # registered factor ceilings give an absolute bound of 2.5 (the raster itself is [0,1])
            weight = float(row['weight'])
            assert 0.0 < weight <= 2.5
            assert row['persistent'] in ('0', '1')
            assert row['corroborated'] in ('0', '1')
            fams[family] += 1
    assert count == c['cloud']['records'] == 349_445
    assert fams == c['cloud']['weighted'] == {'magnetic': 188_878, 'gravity': 160_567}


def test_whole_brief_is_embedded_in_readme_without_omissions():
    readme = (ROOT / 'README.md').read_text()
    embedded = readme.split('<!-- BEGIN USER BRIEF 20261006 -->')[1].split('<!-- END USER BRIEF 20261006 -->')[0].strip()
    assert embedded == (ROOT / 'docs/user-prompt-20261006.md').read_text().strip()
