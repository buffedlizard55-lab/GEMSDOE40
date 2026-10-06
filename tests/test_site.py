from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlparse

DOCS = Path(__file__).resolve().parents[1] / "docs"


class _PageAudit(HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs = []
        self.ids = []

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if "href" in attributes:
            self.hrefs.append(attributes["href"])
        if "id" in attributes:
            self.ids.append(attributes["id"])


def test_static_site_has_no_duplicate_ids_or_broken_local_links():
    pages = sorted(DOCS.rglob("*.html"))
    assert pages
    for page in pages:
        audit = _PageAudit()
        audit.feed(page.read_text(encoding="utf-8"))
        assert len(audit.ids) == len(set(audit.ids)), f"duplicate id in {page}"
        for href in audit.hrefs:
            parsed = urlparse(href)
            if parsed.scheme or href.startswith("//") or not parsed.path:
                continue
            target = (page.parent / parsed.path).resolve()
            assert target.exists(), f"{page}: broken local link {href}"


def test_site_prominently_offers_the_current_download_with_its_real_status():
    """The front page must make the current download obvious and must not imply a validated score."""
    index = (DOCS / "index.html").read_text(encoding="utf-8")
    summary = (DOCS / "executive-summary.html").read_text(encoding="utf-8")
    import json
    manifest = json.loads((DOCS / "data/current-candidate.json").read_text())
    candidate = DOCS / "downloads" / manifest["filename"]
    assert candidate.is_file()
    for text in (index, summary):
        assert f'downloads/{manifest["filename"]}' in text
        assert "not promoted" in text.lower()
        assert "HOLD — DO NOT SUBMIT" in text
    assert "BUILT, AUDITED, NOT PROMOTED" in summary
    assert manifest["slot_eligible"] is False
    assert manifest["organizer_score"] is None
    assert manifest["weekly_submission_used"] is False


def test_site_prominently_marks_retained_artifacts_on_hold_and_offers_their_downloads():
    index = (DOCS / "index.html").read_text(encoding="utf-8")
    summary = (DOCS / "executive-summary.html").read_text(encoding="utf-8")
    import json
    session2 = json.loads((DOCS / "data/session2-artifacts-20261006.json").read_text())
    for name in (session2["filename"], session2["sibling_candidate"]["filename"],
                 session2["h8asa_candidate"]["filename"]):
        assert (DOCS / "downloads" / name).is_file(), name
        assert name in index, name
        assert name in summary, name
    assert "HOLD — DO NOT SUBMIT H40" in index and "HOLD — DO NOT SUBMIT H40" in summary
    assert "HOLD — DO NOT SUBMIT H7 or H40" in index and "HOLD — DO NOT SUBMIT H7 or H40" in summary
    assert "Nothing earlier is cleared either" in index
    assert "Do not upload either version" in summary


def test_h7_and_h40_artifacts_are_offered_beside_an_explicit_hold():
    import hashlib
    import json
    index = (DOCS / "index.html").read_text(encoding="utf-8")
    summary = (DOCS / "executive-summary.html").read_text(encoding="utf-8")
    h7 = "gemsdoe40-h7-rtp-euler-gravity-context-3d-kde-20261006-998f660f.tif"
    assert (DOCS / "downloads" / h7).is_file()
    assert h7 in index and h7 in summary
    zeros = DOCS / "downloads" / "gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-zeros.tif"
    nan = DOCS / "downloads" / "gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-nan.tif"
    receipt_path = DOCS / "downloads" / "gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-audit.json"
    assert zeros.is_file() and nan.is_file() and receipt_path.is_file()
    assert zeros.name in index and zeros.name in summary
    receipt = json.loads(receipt_path.read_text())
    assert hashlib.sha256(zeros.read_bytes()).hexdigest() == receipt["zeros_tif"]["sha256"]
    inst = receipt["stage"]["emission"]["instrument"]
    assert inst["predicted_live"] < 0.2778
    assert inst["lm_calibrated"] < inst["lm_incumbent"]
    assert receipt["stage"]["uniqueness"]["is_new"] is True


def test_site_states_that_no_local_instrument_reaches_the_promotion_bar():
    """The 2026-10-06 audit measured every instrument against 16 owner-attributed scores."""
    import json
    audit = json.loads((DOCS / "data" / "instrument-audit-20261006.json").read_text())
    stats = audit["statistics"]
    assert stats["spearman_lm_vs_live"] <= 0.3        # LM instrument: no ranking power
    assert stats["spearman_mass_vs_live"] <= -0.5     # mass alone anti-correlates
    assert stats["spearman_w_offcat_vs_live"] <= 0.8  # off-catalogue proxy below the bar
    cases = {c["case"] for c in audit["decisive_counterexamples"]}
    assert {"mass-matched control pair", "blind lattice"} <= cases
    assert ">= 0.8" in audit["verdict"] or "0.8" in audit["next_session_requirement"]
    assert len(audit["artifacts"]) == 16
    index = (DOCS / "index.html").read_text(encoding="utf-8")
    evidence = (DOCS / "evidence.html").read_text(encoding="utf-8")
    for text in (index, evidence):
        assert "instrument-audit-20261006.json" in text
    assert "withdrawn" in index.lower()
    assert "≥ 0.800" in evidence
