import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]


def test_real_download_is_float32_and_bounded_with_exact_declared_footprint():
    c = json.loads((ROOT/'docs/data/current-candidate.json').read_text())
    p = ROOT/c['path']
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
    if c["experiment"] == "H13":
        # H13 is the frozen crest-binary emission: exactly {0, 1} values.
        assert np.count_nonzero(a[finite]) == 21_041
        assert np.unique(a[finite]).size == 2
    else:
        raise AssertionError(f"unexpected current experiment {c['experiment']}; pin its measured facts here")


def test_real_cloud_record_count_and_depth_labels_match_measured_generation():
    c = json.loads((ROOT/'docs/data/current-candidate.json').read_text())
    path = ROOT/c['cloud']['path']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == c['cloud']['sha256']
    count, clustered = 0, {'tmi': 0, 'iso_grav_anom': 0}
    with gzip.open(path, 'rt') as stream:
        for row in csv.DictReader(stream):
            count += 1
            family = row['family']
            assert family in clustered
            depth = float(row['effective_depth_m'])
            assert 100 <= depth <= (5000 if family == 'tmi' else 8000)
            assert int(row['identifiable_rank']) in (2,3)
            weight = float(row['cluster_weight'])
            assert 0 <= weight <= 1
            if weight > 0: clustered[family] += 1
    assert count == 46_656
    assert clustered == {'tmi':12_167,'iso_grav_anom':8_450}


def test_whole_brief_is_embedded_in_readme_without_omissions():
    readme=(ROOT/'README.md').read_text()
    embedded=readme.split('<!-- BEGIN USER BRIEF 20261006 -->')[1].split('<!-- END USER BRIEF 20261006 -->')[0].strip()
    assert embedded == (ROOT/'docs/user-prompt-20261006.md').read_text().strip()
