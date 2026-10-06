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


def test_site_prominently_marks_candidate_on_hold_and_offers_download():
    index = (DOCS / "index.html").read_text(encoding="utf-8")
    summary = (DOCS / "executive-summary.html").read_text(encoding="utf-8")
    import json
    manifest = json.loads((DOCS / "data/current-candidate.json").read_text())
    candidate = DOCS / "downloads" / manifest["filename"]
    assert candidate.is_file()
    assert "HOLD — DO NOT SUBMIT" in index
    assert "HOLD — DO NOT SUBMIT" in summary
    assert f'downloads/{manifest["filename"]}' in index
    assert f'downloads/{manifest["filename"]}' in summary
    assert manifest["slot_eligible"] is False
    assert manifest["organizer_score"] is None



def test_site_prominently_marks_h7_on_hold_and_offers_research_download():
    index = (DOCS / "index.html").read_text(encoding="utf-8")
    summary = (DOCS / "executive-summary.html").read_text(encoding="utf-8")
    candidate_name = "gemsdoe40-h7-rtp-euler-gravity-context-3d-kde-20261006-998f660f.tif"
    candidate = DOCS / "downloads" / candidate_name
    assert candidate.is_file()
    assert "HOLD — DO NOT SUBMIT" in index
    assert "HOLD — DO NOT SUBMIT" in summary
    assert candidate_name in index
    assert candidate_name in summary
    assert "Nothing earlier is cleared either" in index
    assert "Do not upload either version" in summary


def test_h40_artifact_is_offered_beside_an_explicit_hold():
    import hashlib
    import json
    index = (DOCS / "index.html").read_text(encoding="utf-8")
    summary = (DOCS / "executive-summary.html").read_text(encoding="utf-8")
    zeros = DOCS / "downloads" / "gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-zeros.tif"
    twin = DOCS / "downloads" / "gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-nan.tif"
    receipt_path = DOCS / "downloads" / "gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-audit.json"
    assert zeros.is_file() and twin.is_file() and receipt_path.is_file()
    assert f"downloads/{zeros.name}" in index and f"downloads/{zeros.name}" in summary
    assert "HOLD — DO NOT SUBMIT" in index
    assert "HOLD — DO NOT SUBMIT H40" in index and "HOLD — DO NOT SUBMIT H40" in summary
    receipt = json.loads(receipt_path.read_text())
    assert hashlib.sha256(zeros.read_bytes()).hexdigest() == receipt["zeros_tif"]["sha256"]
    inst = receipt["stage"]["emission"]["instrument"]
    assert inst["predicted_live"] < 0.2778
    assert inst["lm_calibrated"] < inst["lm_incumbent"]
    assert receipt["stage"]["uniqueness"]["is_new"] is True


def test_instrument_audit_retires_the_two_promotion_instruments():
    """The 2026-10-06 audit measured both instruments against 16 owner-attributed scores; file-level organizer attribution is unverified."""
    import json
    audit = json.loads((DOCS / "data" / "instrument-audit-20261006.json").read_text())
    stats = audit["statistics"]
    assert stats["spearman_lm_vs_live"] <= 0.3        # LM instrument: no ranking power
    assert stats["spearman_mass_vs_live"] <= -0.5     # mass alone anti-correlates
    cases = {c["case"] for c in audit["decisive_counterexamples"]}
    assert {"mass-matched control pair", "blind lattice"} <= cases
    assert ">= 0.8" in audit["verdict"] or "0.8" in audit["next_session_requirement"]
    assert len(audit["artifacts"]) == 16
    index = (DOCS / "index.html").read_text(encoding="utf-8")
    summary = (DOCS / "executive-summary.html").read_text(encoding="utf-8")
    assert "instrument-audit-20261006.json" in index
    assert "instrument-audit-20261006.json" in summary
    assert "withdrawn" in index.lower()
