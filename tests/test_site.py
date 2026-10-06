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
    candidate = DOCS / "downloads/gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-zeros.tif"
    assert candidate.is_file()
    assert "HOLD — DO NOT SUBMIT" in index
    assert "HOLD — DO NOT SUBMIT the 2026-10-06 Euler artifact" in summary
    assert f"downloads/{candidate.name}" in index
    assert f"downloads/{candidate.name}" in summary
    twin = DOCS / "downloads/gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-nan.tif"
    assert twin.is_file()


def test_current_artifact_receipt_shows_no_promotion():
    import hashlib
    import json
    receipt = json.loads((DOCS / "downloads" /
                          "gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-audit.json").read_text())
    zeros = DOCS / "downloads" / "gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-zeros.tif"
    assert hashlib.sha256(zeros.read_bytes()).hexdigest() == receipt["zeros_tif"]["sha256"]
    inst = receipt["stage"]["emission"]["instrument"]
    assert inst["predicted_live"] < 0.2778
    assert inst["lm_calibrated"] < inst["lm_incumbent"]
    assert receipt["stage"]["uniqueness"]["is_new"] is True
