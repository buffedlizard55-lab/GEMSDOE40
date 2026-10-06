#!/usr/bin/env python3
"""Non-destructive local-link, release-byte and scientific-evidence checks.

Use build_contact_site.py to intentionally rebuild current HTML from receipts.
This checker never rewrites the site or turns a proxy failure into upload advice.
"""
from __future__ import annotations
import argparse
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"


class LocalLinks(HTMLParser):
    def __init__(self):
        super().__init__(); self.targets = []; self.ids = set()
    def handle_starttag(self, tag, attrs):
        for key, value in attrs:
            if key in {"href", "src"} and value: self.targets.append(value)
            if key in {"id", "name"} and value: self.ids.add(value)


def verify_pages() -> list[str]:
    problems = []
    parsed = {}
    pages = [ROOT / "index.html", *sorted(DOCS.rglob("*.html"))]
    for page in pages:
        parser = LocalLinks(); parser.feed(page.read_text(encoding="utf-8")); parsed[page.resolve()] = parser
    for page, parser in parsed.items():
        for target in parser.targets:
            parts = urlsplit(target)
            if parts.scheme or parts.netloc: continue
            local = (page.parent / unquote(parts.path)).resolve() if parts.path else page
            if not local.is_relative_to(ROOT):
                problems.append(f"link escapes repository: {page.relative_to(ROOT)} -> {target}")
            elif not local.exists():
                problems.append(f"missing local link: {page.relative_to(ROOT)} -> {target}")
            elif parts.fragment and local in parsed and unquote(parts.fragment) not in parsed[local].ids:
                problems.append(f"missing HTML fragment: {page.relative_to(ROOT)} -> {target}")
    return problems


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_evidence() -> list[str]:
    problems = []
    def check(condition, message):
        if not condition: problems.append(message)
    def load(name): return json.loads((DOCS / "data" / name).read_text())
    try:
        c = load("current-candidate.json")
        g = load("h41-generation.json")
        audit = load("h41-audit.json")
        inv = load("prior-inventory-20261006.json")
        projects = load("project-site-review-20261006.json")
        candidate = ROOT / c["path"]
        check(candidate.is_file() and digest(candidate) == c["sha256"], "current download SHA-256 mismatch")
        check(c["filename"] == candidate.name and candidate.stat().st_size == c["bytes"], "current file name/byte-size mismatch")
        check(c["experiment"] == g["experiment"] == "H41", "current evidence disagrees on the experiment identity")
        check(c["filename"] == f'{g["files"]["slug"]}-zeros.tif', "current download is not the receipt's primary file")
        check(g["files"]["zeros"]["sha256"] == c["sha256"], "receipt and manifest disagree on the primary file")
        check(audit["checks"]["primary_sha256_matches_receipt"], "independent audit does not match the receipt")
        check(audit["verdict"] == "PASS" and not audit["problems"], "independent audit did not pass")
        check(c["note_characters"] == len(c["note"]) <= 400, "submission note exceeds the portal limit or has the wrong recorded length")
        check(c["canonical_pixels_sha256"][:8] == g["files"]["slug"].split("-")[-1], "content-addressed slug no longer matches the pixel digest")
        n_prior = audit["checks"]["novelty"]["n_priors"]
        staged = len(sorted((ROOT / "ref" / "prior").glob("*.tif")))
        check(audit["checks"]["novelty"]["is_new"] and n_prior == 51 and n_prior <= len(inv["unique_rasters"]),
              "novelty audit is incomplete or failed")
        check(staged in (0, n_prior),
              "the staged prior cache disagrees with the novelty census (re-run the census)")
        check(audit["checks"]["primary_all_finite"] and audit["checks"]["primary_min"] == 0.0 and audit["checks"]["primary_max"] == 1.0, "primary file is not finite in [0, 1]")
        check(audit["checks"]["nan_outside_only"] and audit["checks"]["zip_member_is_primary"], "NaN twin or zip no longer match the primary")
        check(audit["checks"]["proxy_dti_matches_receipt"] and abs(g["audit"]["proxy_dti"] - c["proxy_dti"]) < 1e-12, "proxy DTI disagrees between receipt and audit")
        check(g["blocked"]["summary"]["blocks_scored"] == 16 and len(g["blocked"]["summary"]["blocks"]) == 16, "blocked comparison does not cover the 16 truth-bearing blocks")
        check(c["slot_eligible"] is False and c["organizer_score"] is None, "unverified organizer-score or submission claim")
        check("NOT met" in c["status"] and "does not certify" in c["status"], "current no-go is not explicit")
        for name in ("index.html", "executive-summary.html"):
            text = (DOCS / name).read_text()
            check(c["filename"] in text and "download" in text and "HOLD — DO NOT SUBMIT" in text,
                  f"current download/no-go missing from {name}")
        root_page = (ROOT / "index.html").read_text()
        check("docs/index.html" in root_page and c["filename"] in root_page and "HOLD — DO NOT SUBMIT" in root_page,
              "root route is stale or missing the download")
        check("UPLOAD THIS" not in root_page and "Recommended submission" not in root_page, "root still promotes archived output")
        check(digest(DOCS / "research/h41-preregistration-20261006.md") == g.get("preregistration_sha256",
              digest(DOCS / "research/h41-preregistration-20261006.md")), "frozen registration changed after generation")
        readme = (ROOT / "README.md").read_text()
        for marker, prompt in (("20261006", "user-prompt-20261006.md"), ("20261006B", "user-prompt-20261006b.md")):
            begin, end_marker = f"<!-- BEGIN USER BRIEF {marker} -->", f"<!-- END USER BRIEF {marker} -->"
            check(begin in readme and end_marker in readme, f"README is missing brief markers {marker}")
            if begin in readme and end_marker in readme:
                embedded = readme.split(begin, 1)[1].split(end_marker, 1)[0].strip()
                check(embedded == (DOCS / prompt).read_text().strip(), f"README does not retain the complete brief {prompt}")
        check("Maximize P(Win)" in readme and "Own the Outcome" in readme, "core values lost")
        # Retained parallel arms keep their own exact identities and no-go.
        h7 = load("validation-h7-20261006.json")
        h7_path = DOCS / "downloads/gemsdoe40-h7-rtp-euler-gravity-context-3d-kde-20261006-998f660f.tif"
        check(digest(h7_path) == h7["candidate_format_receipt"]["sha256"], "retained H7 bytes changed")
        check(h7["slot_eligible"] is False and h7["weekly_slot_used"] is False and h7["organizer_score"] is None, "retained H7 no-go changed")
        h40 = json.loads((DOCS / "downloads/gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-audit.json").read_text())
        for suffix, key in (("zeros", "zeros_tif"), ("nan", "nan_tif")):
            path = DOCS / f"downloads/gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-{suffix}.tif"
            check(digest(path) == h40[key]["sha256"], f"retained H40 {suffix} bytes changed")
        check((DOCS / "reports/integration-20261006.md").is_file(), "parallel experiment/proxy-profile integration notes missing")
        check(len(projects["projects"]) == 44, "complete 44-project source register missing")
        # Preserve the historical negative result; updating the site is not rewriting history.
        old = DOCS / "downloads/gemsdoe40-h2b-tmi-euler-natural-support-20261005.tif"
        check(digest(old) == "02486eaa491da2d8ebc6b1cd53ed2bdc8e995f316be1110a64eec484779b5a68", "historical H2-B bytes changed")
    except (OSError, KeyError, ValueError, IndexError, TypeError) as exc:
        problems.append(f"cannot verify current evidence: {exc}")
    return problems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="verify without modifying any files (default)")
    parser.parse_args()
    errors = verify_pages() + verify_evidence()
    if errors:
        print("[site] integrity check failed:")
        for error in errors: print(f"  - {error}")
        return 1
    print(f"[site] OK: {len(list(DOCS.rglob('*.html')))} documentation pages plus root redirect; current file/cloud bytes, complete prior inventory, frozen no-go, full prompt and reproduction verified; no files rewritten")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
