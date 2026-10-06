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


def test_site_prominently_marks_h7_on_hold_and_offers_research_download():
    index = (DOCS / "index.html").read_text(encoding="utf-8")
    summary = (DOCS / "executive-summary.html").read_text(encoding="utf-8")
    candidate_name = "gemsdoe40-h7-rtp-euler-gravity-context-3d-kde-20261006-998f660f.tif"
    candidate = DOCS / "downloads" / candidate_name
    assert candidate.is_file()
    assert "HOLD — DO NOT SUBMIT" in index
    assert "HOLD — DO NOT SUBMIT H7" in summary
    assert candidate_name in index
    assert candidate_name in summary
    assert "No candidate is cleared today" in summary
    assert "Do not upload H7" in summary
