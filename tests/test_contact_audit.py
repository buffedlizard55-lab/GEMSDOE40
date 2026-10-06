import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

from gemsdoe40.contact_audit import audit, compare_vectors, correlation, top_indices
from gemsdoe40.raster import write_candidate


def test_ties_are_reproducible_and_zeros_not_selected():
    np.testing.assert_array_equal(top_indices(np.array([0, 1, 1, 1, 0]), 2), [1, 2])
    assert not len(top_indices(np.zeros(20), 4))


def test_rescaling_or_tiny_noise_is_not_novel():
    x = np.linspace(0, 1, 1000, dtype=np.float32)
    for y in (x.copy(), x*.7, x + 1e-5*np.sin(np.arange(x.size))):
        assert compare_vectors(x, y, budget=100)["near_duplicate_h4"]


def test_different_continuous_maps_do_not_fail_merely_for_shared_background():
    rng = np.random.default_rng(4)
    x, y = rng.uniform(.001, 1, (2, 10000)).astype(np.float32)
    r = compare_vectors(x, y, budget=100)
    assert r["nonzero_support_jaccard"] == 1
    assert r["historic_h2b_rule_near_duplicate"]
    assert not r["near_duplicate_h4"]


def test_constant_correlation_and_invalid_vectors():
    assert correlation(np.ones(5), np.zeros(5)) == 0
    assert correlation(np.zeros(5), np.zeros(5)) == 1
    with pytest.raises(ValueError, match="nonempty"):
        compare_vectors(np.array([]), np.array([]))
    with pytest.raises(ValueError, match="finite"):
        compare_vectors(np.array([np.nan]), np.array([.2]))


def files(tmp_path):
    sample, candidate = tmp_path / "sample.tif", tmp_path / "candidate.tif"
    s = np.zeros((12, 12), np.float32); s[-1] = np.nan
    with rasterio.open(sample, "w", driver="GTiff", count=1, dtype="float32", nodata=np.nan,
                       width=12, height=12, crs="EPSG:32611", transform=from_origin(243350, 4508550, 100, 100)) as ds:
        ds.write(s, 1)
    p = np.linspace(0, 1, 144, dtype=np.float32).reshape(12, 12)
    write_candidate(candidate, p, sample)
    cache = tmp_path / "cache"; cache.mkdir()
    return sample, candidate, cache


def test_missing_prior_is_fail_closed(tmp_path):
    sample, candidate, cache = files(tmp_path)
    inventory = tmp_path / "inv.json"
    inventory.write_text(json.dumps({"unique_rasters": [{"git_blob_sha": "absent", "bytes": 1, "sha256": "absent"}]}))
    report = audit(candidate, sample, inventory, cache)
    assert not report["uniqueness_pass"]
    assert not report["completeness_pass"]
    assert report["issues"]


def test_empty_inventory_is_not_a_uniqueness_pass(tmp_path):
    sample, candidate, cache = files(tmp_path)
    inventory = tmp_path / "inv.json"; inventory.write_text('{"unique_rasters": []}')
    with pytest.raises(ValueError, match="empty/incomplete"):
        audit(candidate, sample, inventory, cache)


def test_exact_raw_copy_veto_and_integrity_check(tmp_path):
    sample, candidate, cache = files(tmp_path)
    data = candidate.read_bytes()
    blob = hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()
    p = cache / f"{blob}.tif"; p.write_bytes(data)
    entry = {"git_blob_sha": blob, "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data),
             "artifacts": [{"path": "candidate.tif", "repo": "fixture", "ref": "test"}]}
    inventory = tmp_path / "inv.json"; inventory.write_text(json.dumps({"unique_rasters": [entry]}))
    result = audit(candidate, sample, inventory, cache)
    assert result["completeness_pass"]
    assert not result["uniqueness_pass"]
    assert result["comparisons"][0]["file_hash_equal"]
    assert result["comparisons"][0]["canonical_hash_equal"]
    p.write_bytes(data + b"corrupt")
    with pytest.raises(ValueError, match="corrupt prior"):
        audit(candidate, sample, inventory, cache)


def _entry_bytes(data, cache, name="prior.tif"):
    blob = hashlib.sha1(f"blob {len(data)}\0".encode()+data).hexdigest()
    (cache/f"{blob}.tif").write_bytes(data)
    return {"git_blob_sha": blob, "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data),
            "artifacts": [{"path": name, "repo": "fixture", "ref": "test"}]}


def test_unclassified_unreadable_prior_still_fails_closed(tmp_path):
    sample, candidate, cache = files(tmp_path)
    entry = _entry_bytes(b"not actually a TIFF, even with an honest hash", cache)
    inventory = tmp_path/"inv.json"; inventory.write_text(json.dumps({"unique_rasters": [entry]}))
    result = audit(candidate, sample, inventory, cache)
    assert not result["completeness_pass"] and not result["uniqueness_pass"]
    assert "unclassified unreadable" in result["issues"][0]["issue"]


def test_exact_header_only_artifact_is_accounted_not_falsely_correlated(tmp_path):
    sample, candidate, cache = files(tmp_path)
    # Actual measured 110-byte historical artifact: IFD ends at EOF; no pixels.
    data = bytes.fromhex("49492a000800000008000001030001000000dc0c00000101030001000000920e00000201030001000000200000000301030001000000050000000601030001000000010000001501030001000000010000001c010300010000000200000053010300010000000300000000000000")
    header = _entry_bytes(data, cache)
    inventory = tmp_path/"inv.json"; inventory.write_text(json.dumps({"unique_rasters": [header]}))
    result = audit(candidate, sample, inventory, cache)
    assert result["nonprediction_artifacts_forensically_classified"] == 1
    assert result["hashed_and_audited_blobs"] == 1
    assert not result["comparisons"][0]["correlation_applicable"]
    assert "pearson_r" not in result["comparisons"][0]
    # A corpus of only technical fixtures is not a valid novelty claim.
    assert not result["uniqueness_pass"] and not result["completeness_pass"]


def test_fixture_filename_is_not_a_blanket_novelty_exemption(tmp_path):
    sample, candidate, cache = files(tmp_path)
    file = tmp_path/"small.tif"
    with rasterio.open(file, "w", driver="GTiff", count=1, dtype="float32", width=2, height=2,
                       crs="EPSG:32611", transform=from_origin(243350, 4508550, 100, 100)) as ds:
        ds.write(np.zeros((2,2), np.float32),1)
    entry = _entry_bytes(file.read_bytes(), cache, "FORMAT_TEST_NOT_SUBMISSION.tif")
    inventory = tmp_path/"inv.json"; inventory.write_text(json.dumps({"unique_rasters": [entry]}))
    result = audit(candidate, sample, inventory, cache)
    assert not result["uniqueness_pass"]
    assert any("exact-byte known" in i["issue"] for i in result["issues"])


def test_nan_storage_equivalence_requires_exact_real_companion_values():
    from gemsdoe40.contact_audit import nan_zero_storage_equivalent
    partial=np.array([[.4,np.nan],[0,.8]],np.float32)
    complete=np.array([[.4,0],[0,.8]],np.float32)
    assert nan_zero_storage_equivalent(partial,complete)
    changed=complete.copy();changed[0,1]=.001
    assert not nan_zero_storage_equivalent(partial,changed)
    assert not nan_zero_storage_equivalent(partial,partial)
    assert not nan_zero_storage_equivalent(np.array([np.inf]),np.array([0.]))


def test_low_coverage_variant_requires_actual_pinned_companion_in_complete_corpus(tmp_path,monkeypatch):
    import gemsdoe40.contact_audit as module
    sample,candidate,cache=files(tmp_path)
    def raw_file(name,values):
        path=tmp_path/name
        with rasterio.open(path,"w",driver="GTiff",count=1,dtype="float32",nodata=np.nan,
                           width=12,height=12,crs="EPSG:32611",transform=from_origin(243350,4508550,100,100)) as ds:
            ds.write(values,1)
        return _entry_bytes(path.read_bytes(),cache,name)
    full=raw_file('full.tif',np.zeros((12,12),np.float32))
    a=np.zeros((12,12),np.float32);a[2:]=np.nan
    partial=raw_file('partial.tif',a)
    inventory=tmp_path/'inv.json';inventory.write_text(json.dumps({'unique_rasters':[partial,full]}))
    monkeypatch.setattr(module,'NAN_STORAGE_COUNTERPARTS',{partial['sha256']:(full['git_blob_sha'],full['sha256'])})
    result=audit(candidate,sample,inventory,cache)
    assert result['completeness_pass'] and result['uniqueness_pass']
    row=result['comparisons'][0]
    assert row['comparison_coverage']<.9
    assert row['pixels_compared']==24  # Do not silently impute the missing pixels.
    assert row['verified_storage_counterpart']['counterpart_nonzero_at_missing_cells']==0
    (cache/f"{full['git_blob_sha']}.tif").unlink()
    result=audit(candidate,sample,inventory,cache)
    assert not result['uniqueness_pass'] and not result['completeness_pass']
