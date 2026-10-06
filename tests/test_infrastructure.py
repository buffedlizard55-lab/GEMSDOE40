from __future__ import annotations
from datetime import datetime, timezone
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import zipfile

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


def script(name):
    spec = importlib.util.spec_from_file_location(f"tested_{name}", ROOT / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_tif(path, arr, transform=None):
    with rasterio.open(path, "w", driver="GTiff", width=arr.shape[1], height=arr.shape[0], count=1,
                       dtype="float32", crs="EPSG:32611", transform=transform or from_origin(243350, 4508550, 100, 100), nodata=np.nan) as ds:
        ds.write(arr.astype(np.float32), 1)


def test_prepare_refuses_missing_inputs_without_overwriting_pins(tmp_path, monkeypatch):
    module = script("prepare_data"); monkeypatch.setattr(module, "ROOT", tmp_path)
    with pytest.raises(SystemExit, match="Missing required inputs"):
        module.main()
    assert not (tmp_path / "evidence/data_pins.json").exists()


def test_prepare_refuses_shifted_geotransform(tmp_path, monkeypatch):
    module = script("prepare_data"); monkeypatch.setattr(module, "ROOT", tmp_path)
    data = tmp_path / "data"; data.mkdir()
    a = np.zeros((4, 5), dtype=np.float32)
    for name in module.PINS:
        write_tif(data / name, a, from_origin(243450, 4508550, 100, 100) if name == "labels.tif" else None)
    monkeypatch.setattr(module, "PINS", {name: module.sha256(data/name) for name in module.PINS})
    with pytest.raises(SystemExit, match="unaligned input"):
        module.main()
    assert not (tmp_path / "evidence/data_pins.json").exists()


def test_acquisition_is_atomic_and_checks_bytes(tmp_path, monkeypatch):
    module = script("acquire_data")
    payload = b"a fixed test payload, not a raster"
    calls = []
    def fake_run(command, **kwargs):
        calls.append(command); kwargs["stdout"].write(payload)
    monkeypatch.setattr(module.subprocess, "run", fake_run)
    dest = tmp_path / "source.tif"
    module.fetch("test", dest, sha256=hashlib.sha256(payload).hexdigest(), size=len(payload))
    module.fetch("test", dest, sha256=hashlib.sha256(payload).hexdigest(), size=len(payload))
    assert dest.read_bytes() == payload and len(calls) == 1
    assert not list(tmp_path.glob("*.partial"))
    wrong = tmp_path / "wrong.tif"
    with pytest.raises(RuntimeError, match="integrity check"):
        module.fetch("test", wrong, sha256="0"*64)
    assert not wrong.exists() and not list(tmp_path.glob("*.partial"))


def test_acquisition_preserves_an_unexpected_existing_file(tmp_path):
    module = script("acquire_data")
    dest = tmp_path / "source.tif"; dest.write_bytes(b"do not overwrite")
    with pytest.raises(RuntimeError, match="preserving"):
        module.fetch("not-requested", dest, sha256="0"*64)
    assert dest.read_bytes() == b"do not overwrite"


def test_prior_zip_member_read_is_not_filesystem_extraction(tmp_path):
    module = script("fetch_prior_cache")
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("../never-extract.tif", b"verified TIFF bytes")
    result = module.unwrap_artifact(buffer.getvalue(), {"member": "../never-extract.tif"})
    assert result == b"verified TIFF bytes"
    assert not (tmp_path.parent / "never-extract.tif").exists()
    with pytest.raises(KeyError):
        module.unwrap_artifact(buffer.getvalue(), {"member": "missing.tif"})


def test_prior_raw_payload_is_not_changed():
    module = script("fetch_prior_cache")
    assert module.unwrap_artifact(b"raw", {}) == b"raw"
    assert module.git_blob_sha1(b"raw") == hashlib.sha1(b"blob 3\0raw").hexdigest()


def test_classifier_recovers_raw_and_historical_prediction_locations():
    module = script("refresh_prior_inventory")
    for path in ("data/prob_raw.tif", "data/context_detector_prob.tif", "inputs/calibration/example.tif",
                 "data/evidence/leaderboard_anchor/old.tif", "data/evidence/suture/ens_only.tif",
                 "data/evidence/union_po_loo/prev_committed_nff.tif", "inputs/gems19-old-nan.tif",
                 "docs/archive/legacy.zip", "docs/gems40-output.tif"):
        assert module.output_reason(path), path
    for path in ("data/labels.tif", "data/sample_submission.tif", "external/qfaults_prior_u8.tif"):
        assert module.output_reason(path) is None


def test_every_prompt_score_is_retained_in_the_register():
    module = script("review_project_sites")
    parsed = module.parse_prompt((ROOT / "docs/user-prompt-20261006.md").read_text())
    register = json.loads((ROOT / "docs/data/project-site-review-20261006.json").read_text())
    assert len(parsed) == 44
    assert sum(s["score"] is not None for p in parsed for s in p["submissions"]) == 48
    assert [(p["repo"], p["submissions"]) for p in parsed] == [(p["repo"], p["submissions"]) for p in register["projects"]]
    assert any(not p["submissions"] for p in parsed)


def test_context_feed_is_not_a_drivendata_scraper():
    module = script("update_source_feed")
    assert all("drivendata" not in host for host in module.ALLOWED_HOSTS)
    url = module.earthquake_url(datetime(2026, 10, 6, tzinfo=timezone.utc))
    assert "earthquake.usgs.gov" in url and "2026-09-06" in url and "limit=20" in url
    with pytest.raises(ValueError, match="restricted"):
        module.get("https://www.drivendata.org/competitions/306/leaderboard/")


def test_feed_preserves_last_good_context_and_marks_network_failure(tmp_path, monkeypatch):
    module = script("update_source_feed")
    out = tmp_path / "feed.json"
    old = {"last_success_utc": "2026-10-05T10:00:00+00:00", "events": [{"id": "old-context"}]}
    out.write_text(json.dumps(old))
    def unavailable(url): raise OSError("TLS/network unavailable")
    monkeypatch.setattr(module, "get", unavailable)
    result = module.update(out)
    assert result["status"] == "stale" and result["last_success_utc"] == old["last_success_utc"]
    assert result["events"] == old["events"] and result["errors"]
    assert result["source_probe_failures"] == 3
    assert not out.with_suffix(".tmp").exists()


def test_feed_valid_empty_result_is_not_fabricated_failure(tmp_path, monkeypatch):
    module = script("update_source_feed")
    monkeypatch.setattr(module, "get", lambda url: (b'{"type":"FeatureCollection","features":[]}', 200))
    result = module.update(tmp_path / "feed.json")
    assert result["status"] == "ok" and result["events"] == [] and result["last_success_utc"]
    assert not result["errors"] and result["source_probe_failures"] == 0


def test_feed_rejects_non_usgs_event_link():
    module = script("update_source_feed")
    data = {"type": "FeatureCollection", "features": [{"id": "x", "properties": {"url": "https://example.org/"}}]}
    with pytest.raises(ValueError, match="event-link host"):
        module.parse_events(json.dumps(data).encode())


def test_legacy_runners_require_deliberate_opt_in(monkeypatch, capsys):
    from gemsdoe40.legacy_guard import require_legacy_opt_in
    monkeypatch.delenv("GEMSDOE_ALLOW_LEGACY_RESEARCH", raising=False)
    with pytest.raises(SystemExit, match="ARCHIVED RESEARCH"):
        require_legacy_opt_in()
    monkeypatch.setenv("GEMSDOE_ALLOW_LEGACY_RESEARCH", "1")
    require_legacy_opt_in()
    assert "not approval to submit" in capsys.readouterr().err


def test_legacy_empty_directory_is_not_uniqueness(tmp_path):
    from gemsdoe40.uniqueness import compare_against, pearson
    report = compare_against(np.zeros((4, 5)), tmp_path)
    assert not report["is_new"] and not report["comparison_complete"]
    assert pearson(np.zeros((4, 5)), np.zeros((4, 5))) == 1


def test_legacy_inventory_audit_fails_closed_on_missing_cache(tmp_path):
    from gemsdoe40.research_uniqueness import audit_candidate
    a = np.zeros((4, 5), dtype=np.float32)
    for name in ("sample.tif", "candidate.tif"): write_tif(tmp_path/name, a)
    inv = tmp_path / "inventory.json"; inv.write_text(json.dumps({"unique_rasters": [{"git_blob_sha": "f"*40, "sha256": "0"*64}]}))
    result = audit_candidate(tmp_path/"candidate.tif", tmp_path/"sample.tif", inv, tmp_path/"absent")
    assert not result["uniqueness_pass"] and not result["completeness_pass"]


def test_generation_cannot_write_into_published_docs():
    result = subprocess.run([sys.executable, str(ROOT/"scripts/run_contact_euler.py"), "--output", str(ROOT/"docs/downloads")], capture_output=True, text=True)
    assert result.returncode != 0 and "ignored work/" in result.stderr
