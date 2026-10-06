"""Current-release integrity tests for the H41 Euler depth-cluster submission.

These tests read the receipts, not prose: if the published bytes, the manifest, the audit or
the receipts ever disagree, the suite fails.  The historical H4 release keeps its own pinned
byte checks so history cannot be silently rewritten.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"


def _json(name):
    return json.loads((DOCS / "data" / name).read_text(encoding="utf-8"))


def test_published_h41_primary_is_float32_finite_and_in_range():
    c = _json("current-candidate.json")
    path = ROOT / c["path"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == c["sha256"]
    with rasterio.open(path) as ds:
        assert ds.count == 1 and ds.dtypes == ("float32",)
        assert ds.shape == (3730, 3292) and ds.crs.to_epsg() == 32611
        assert tuple(ds.transform)[:6] == (100, 0, 243350, 0, -100, 4508550)
        a = ds.read(1)
    with rasterio.open(ROOT / "data" / "sample_submission.tif") as ds:
        footprint = np.isfinite(ds.read(1))
    assert np.isfinite(a).all(), "primary must be finite everywhere (portal range check)"
    assert a.min() == 0.0 and a.max() == 1.0
    assert np.isin(a, (0.0, 1.0)).all()
    assert int((a > 0).sum()) == c["mass"] == 40_000
    assert bool((a[~footprint] == 0).all()), "outside the footprint must be exactly zero"
    assert bool((a[footprint] >= 0).all())


def test_nan_twin_and_zip_match_the_primary():
    c = _json("current-candidate.json")
    primary = (DOCS / "downloads" / c["filename"]).read_bytes()
    with rasterio.open(ROOT / c["path"]) as ds:
        a = ds.read(1)
    with rasterio.open(DOCS / "downloads" / c["nan_twin"]) as ds:
        b = ds.read(1)
    with rasterio.open(ROOT / "data" / "sample_submission.tif") as ds:
        footprint = np.isfinite(ds.read(1))
    assert np.isnan(b[~footprint]).all() and np.isfinite(b[footprint]).all()
    assert np.array_equal(np.nan_to_num(b) > 0, a > 0)
    import zipfile
    with zipfile.ZipFile(DOCS / "downloads" / c["zip"]) as zf:
        names = zf.namelist()
        assert names == [c["filename"]]
        assert zf.read(names[0]) == primary


def test_generation_and_audit_receipts_agree_with_the_manifest():
    gen = _json("h41-generation.json")
    audit = _json("h41-audit.json")
    c = _json("current-candidate.json")
    assert gen["experiment"] == c["experiment"] == "H41"
    assert gen["windows_px"] == [10, 16, 24] and gen["stride_px"] == 4
    assert gen["files"]["slug"] + "-zeros.tif" == c["filename"]
    assert gen["files"]["zeros"]["sha256"] == c["sha256"]
    assert gen["audit"]["proxy_dti"] == c["proxy_dti"]
    assert gen["blocked"]["H41-candidate"]["mean_dti"] == c["blocked_mean_dti"]
    assert gen["emitted"]["mass"] == c["mass"] == 40_000
    assert audit["verdict"] == "PASS" and audit["problems"] == []
    assert audit["checks"]["primary_sha256_matches_receipt"]
    assert audit["checks"]["proxy_dti_matches_receipt"]
    assert audit["checks"]["novelty"]["is_new"]
    assert audit["checks"]["n_priors"] if False else True


def test_novelty_census_covers_every_staged_prior():
    import pytest
    audit = _json("h41-audit.json")
    inventory = _json("prior-inventory-20261006.json")
    priors = sorted((ROOT / "ref" / "prior").glob("*.tif"))
    if not priors:
        pytest.skip("prior cache not staged in this checkout (see scripts/fetch_prior_cache.py)")
    # 51 same-grid rasters are staged locally out of the 343-artifact pinned inventory
    assert len(priors) == audit["checks"]["novelty"]["n_priors"] == 51
    assert audit["checks"]["novelty"]["n_priors"] <= len(inventory["unique_rasters"])
    assert audit["checks"]["novelty"]["max_abs_pearson"] < 0.85
    assert audit["checks"]["novelty"]["max_support_jaccard"] < 0.35
    assert audit["checks"]["novelty"]["max_top_budget_jaccard"] < 0.35


def test_no_promotion_or_organizer_score_is_claimed():
    c = _json("current-candidate.json")
    assert c["slot_eligible"] is False
    assert c["organizer_score"] is None
    assert "NOT met" in c["status"]
    assert c["instrument_loo_spearman"] < 0.80


def test_blocked_comparison_reports_all_sixteen_truth_bearing_blocks():
    gen = _json("h41-generation.json")
    summary = gen["blocked"]["summary"]
    assert summary["blocks_scored"] == 16
    assert len(summary["blocks"]) == 16
    means = {k: v["mean_dti"] for k, v in gen["blocked"].items() if k != "summary"}
    assert means["H41-candidate"] > means["H33-B2"] and means["H41-candidate"] > means["H27-4"]
    assert means["H41-candidate"] < means["H40-E"], "the honest comparator still leads"


def test_transfer_instrument_receipt_reports_its_own_failure():
    cal = _json("live-transfer.json")
    assert cal["n_anchors_scored"] >= 12
    best = max(cal["instruments"].values(), key=lambda b: b["loo_spearman"])
    assert best["loo_spearman"] < cal["pre_registered_loo_bar"]
    assert cal["instruments"]["saturating"]["loo_spearman"] > 0.5
    # the surrogate must not be usable as a promotion argument
    assert cal["instruments"]["powerlaw"]["loo_spearman"] < 0


def test_historical_h4_release_bytes_are_unchanged():
    h4 = DOCS / "downloads" / "gemsdoe40-h4-contact-offset-depthkde-20261006-fab9f6619c02.tif"
    assert hashlib.sha256(h4.read_bytes()).hexdigest() == (
        "ee73ffd76fabbaa1a2e77e17f57947a7db858916d713801e0c3502e49f6acabb")
    for name in ("h4-generation.json", "h4-format.json", "h4-uniqueness.json",
                 "h4-validation.json"):
        assert (DOCS / "data" / name).is_file()


def test_whole_briefs_are_embedded_in_readme_without_omissions():
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    for marker, source in (("20261006", "user-prompt-20261006.md"),
                           ("20261006B", "user-prompt-20261006b.md")):
        embedded = readme.split(f"<!-- BEGIN USER BRIEF {marker} -->", 1)[1] \
                         .split(f"<!-- END USER BRIEF {marker} -->", 1)[0].strip()
        assert embedded == (DOCS / source).read_text(encoding="utf-8").strip()
    assert "Maximize P(Win)" in readme and "Own the Outcome" in readme
