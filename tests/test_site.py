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
    import json
    index = (DOCS / "index.html").read_text(encoding="utf-8")
    summary = (DOCS / "executive-summary.html").read_text(encoding="utf-8")
    manifest = json.loads((DOCS / "data/current-candidate.json").read_text())
    candidate = DOCS / "downloads" / manifest["filename"]
    assert candidate.is_file()
    assert "HOLD — DO NOT SUBMIT" in index
    assert "HOLD — DO NOT SUBMIT" in summary
    assert f'downloads/{manifest["filename"]}' in index
    assert f'downloads/{manifest["filename"]}' in summary
    assert manifest["slot_eligible"] is False
    assert manifest["organizer_score"] is None
    # the gate reason must be stated on the pages, not only in the manifest
    assert "0.80" in index or "0.80" in summary
    assert "H40-E" in index or "H40-E" in summary


def test_site_states_the_five_verified_representation_and_submission_facts():
    index = (DOCS / "index.html").read_text(encoding="utf-8")
    summary = (DOCS / "executive-summary.html").read_text(encoding="utf-8")
    evidence = (DOCS / "evidence.html").read_text(encoding="utf-8")
    for text in (index, summary, evidence):
        assert "EPSG:32611" in text
        assert "float32" in text
    assert "356,650" in index            # Euler solutions
    assert "SI = 0" in index or "SI0" in index
    assert "[0, 1]" in summary or "range test" in summary
    assert "Reid" in (DOCS / "sources.html").read_text(encoding="utf-8")


def test_hypotheses_page_lists_five_registered_hypotheses_with_external_dependencies():
    page = (DOCS / "hypotheses.html").read_text(encoding="utf-8")
    for code in ("H41", "H42", "H43", "H44", "H45"):
        assert code in page
    assert "Blocked by" in page            # the GeoDAWN / DEM dependency is disclosed
    assert "surface" in page.lower()


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



def test_instrument_audits_are_retained_and_the_current_one_reports_its_failure():
    """The 2026-10-06 audit retired that session's instruments; the H41 transfer receipt
    publishes the current four models and their leave-one-out failures."""
    import json
    audit = json.loads((DOCS / "data" / "instrument-audit-20261006.json").read_text())
    stats = audit["statistics"]
    assert stats["spearman_lm_vs_live"] <= 0.3        # LM instrument: no ranking power
    assert stats["spearman_mass_vs_live"] <= -0.5     # mass alone anti-correlates
    cases = {c["case"] for c in audit["decisive_counterexamples"]}
    assert {"mass-matched control pair", "blind lattice"} <= cases
    assert len(audit["artifacts"]) == 16
    transfer = json.loads((DOCS / "data" / "live-transfer.json").read_text())
    assert transfer["pre_registered_loo_bar"] == 0.8
    assert all(not blob["passes_bar"] for blob in transfer["instruments"].values())
    evidence = (DOCS / "evidence.html").read_text(encoding="utf-8")
    assert "instrument-audit-20261006.json" in evidence
    assert "live-transfer.json" in evidence
    assert "withdrawn" in evidence.lower()

def test_static_contract_matches_the_browser_qa_expectations():
    """Encode the browser QA's static requirements so they fail in pytest, not only in Chromium.

    scripts/check_site_browser.py cannot run in this sandbox (cdn.playwright.dev is unreachable),
    so the same assertions are checked here against the generated HTML.
    """
    import json
    import re
    release = json.loads((DOCS / "data/current-candidate.json").read_text())
    pages = ("index.html", "executive-summary.html", "evidence.html", "hypotheses.html",
             "methodology.html", "leaderboard-analysis.html", "leaderboard.html", "sources.html")
    for name in pages:
        page = (DOCS / name).read_text(encoding="utf-8")
        assert len(re.findall(r"<h1[ >]", page)) == 1, f"{name}: exactly one h1 required"
        assert re.search(r'<a[^>]*href="#main"[^>]*class="skip"|<a[^>]*class="skip"[^>]*href="#main"', page), \
            f"{name}: skip link"
    for name in ("index.html", "executive-summary.html"):
        page = (DOCS / name).read_text(encoding="utf-8")
        tags = [t for t in re.findall(r"<a[^>]*>", page)
                if f'href="downloads/{release["filename"]}"' in t]
        assert tags and all("download" in t for t in tags), f"{name}: download anchor"
    summary = (DOCS / "executive-summary.html").read_text(encoding="utf-8")
    for control in ("filename", "submission-name", "submission-note", "file-sha"):
        assert f'data-copy="{control}"' in summary and f'id="{control}"' in summary, control
    board = (DOCS / "leaderboard.html").read_text(encoding="utf-8")
    assert 'id="project-search"' in board
    assert board.count("data-project-row") == 44, "the project register must expose 44 rows"
    index = (DOCS / "index.html").read_text(encoding="utf-8")
    marker = index.index(f'downloads/{release["filename"]}')
    assert marker < index.index("<h2"), "the download must precede the first section heading (first viewport)"
    assert marker < 2500, "the download link must sit in the hero, not below it"

WRAP_ANYWHERE = {"mono", "micro", "file-details", "notice", "lead", "prose", "button", "codebox",
                 "tablewrap", "card", "callout", "stat", "eventlist", "submissions", "sourcecard",
                 "label", "tag", "textlink", "formula", "sectionhead", "hero", "kicker", "footer"}
SCROLLS = {"tablewrap", "codebox", "formula", "code"}


class _TokenAudit(HTMLParser):
    """Collect long unbreakable text tokens with the class chain that contains them."""

    def __init__(self):
        super().__init__()
        self.stack = []
        self.risky = []

    def handle_starttag(self, tag, attrs):
        classes = set((dict(attrs).get("class") or "").split())
        self.stack.append(classes)

    def handle_endtag(self, tag):
        if self.stack:
            self.stack.pop()

    def handle_data(self, data):
        for token in data.split():
            if len(token) < 34:
                continue
            chain = set().union(*self.stack) if self.stack else set()
            if chain & WRAP_ANYWHERE or chain & SCROLLS:
                continue
            self.risky.append((token[:60], sorted(chain)))


def test_no_long_unbreakable_token_can_widen_a_page_at_320px():
    """The CI browser leg failed at 320 px once already; keep the cause out of the pages.

    Long single tokens (a 64-char hash, a 37-char filename) only stay inside a 320 px
    viewport when an ancestor wraps anywhere or scrolls horizontally.
    """
    offenders = []
    for page in [DOCS / name for name in (
            "index.html", "executive-summary.html", "evidence.html", "hypotheses.html",
            "methodology.html", "leaderboard-analysis.html", "leaderboard.html", "sources.html")]:
        audit = _TokenAudit()
        audit.feed(page.read_text(encoding="utf-8"))
        offenders.extend((page.name, token, classes) for token, classes in audit.risky)
    assert not offenders, f"unbreakable tokens outside wrapping contexts: {offenders[:5]}"
    # The 2026-10-06 320 px failure was intrinsic sizing of the numeric checkstrip, not a token
    # of this length, so pin the two rules that shrink it.
    css = (DOCS / "assets/contact.css").read_text(encoding="utf-8")
    assert ".checkstrip>div,.checkstrip>*{min-width:0}" in css
    assert ".checkstrip .stat{font-size:clamp(20px,7.4vw,31px);overflow-wrap:anywhere}" in css
