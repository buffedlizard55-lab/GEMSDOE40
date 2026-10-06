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
    """Verify the *current* release chain: H8 candidate, audit and retained history."""
    problems: list[str] = []
    def check(condition, message):
        if not condition:
            problems.append(message)
    def load(name):
        return json.loads((DOCS / "data" / name).read_text())
    try:
        c = load("current-candidate.json")
        audit = load("h8-lineament-audit.json")
        generation = load("h8-lineament-generation.json")
        inv = load("prior-inventory-20261006.json")
        projects = load("project-site-review-20261006.json")

        candidate = ROOT / c["path"]
        check(candidate.is_file(), f"current download missing: {c['path']}")
        check(digest(candidate) == c["sha256"], "current download SHA-256 mismatch")
        check(candidate.name == c["filename"] and candidate.stat().st_size == c["bytes"],
              "current file name/byte-size mismatch")
        check(c["canonical_pixels_sha256"] == audit["candidate"]["canonical_pixels_sha256"],
              "canonical pixel identity disagrees with the audit receipt")
        check(audit["candidate"]["valid"] is True, "current file is not format-verified")
        check(c["sha256"] == audit["published"]["sha256"] == audit["candidate"]["sha256"],
              "publication and audit disagree on the file hash")
        check(audit["novelty"]["novel"] is True and not audit["novelty"]["exact_duplicates"],
              "novelty gate did not pass")
        check(audit["novelty"]["max_abs_pearson"] <= audit["novelty"]["thresholds"]["max_abs_pearson"],
              "raw correlation exceeds the registered novelty threshold")
        check(audit["emitted_dots"] > 0 and audit["mass"] > 0, "empty emission")
        check(audit["emitted_dots"] == c["dots"] == audit["candidate"]["candidate"]["in_footprint_nonzero"],
              "emitted dot count disagrees between receipt and audit")
        check(abs(c["mass"] - audit["mass"]) < 1e-6, "emitted mass disagrees between receipt and audit")
        check(generation["construction_reads_labels_or_priors"] is False,
              "generation claims to read labels/priors")
        check(generation["configuration"]["si"] == 0,
              "registered structural index changed after generation")
        check(abs(generation["configuration"]["nms_spacing_px"] - 2.8) < 1e-9,
              "registered emission spacing changed after generation")
        cloud = ROOT / c["cloud"]["path"]
        check(cloud.is_file() and digest(cloud) == c["cloud"]["sha256"],
              "downloadable depth-labelled solution cloud changed")
        check(c["note_characters"] == len(c["note"]) <= 200, "submission note length is wrong")
        check(c["slot_eligible"] is False and c["organizer_score"] is None
              and c["weekly_submission_used"] is False, "release status became an unverified claim")

        cached = len(list((ROOT / "data" / "prior").glob("*.tif")))
        comparable = audit["novelty"]["corpus_size"] + audit["novelty"]["unreadable_or_incomparable"]
        check(inv["unique_git_blobs"] == len(inv["unique_rasters"]) >= cached,
              "prior inventory is smaller than the cached corpus")
        check(cached == comparable, "cached prior corpus disagrees with the novelty audit")
        check(projects["project_count"] == len(projects["projects"]) == 44,
              "complete 44-project source register missing")
        check(projects["reported_score_count"]
              == sum(s["score"] is not None for p in projects["projects"] for s in p["submissions"]) == 48,
              "48 owner-reported scores were not preserved")

        index_text = (DOCS / "index.html").read_text()
        summary_text = (DOCS / "executive-summary.html").read_text()
        for name, text in (("index.html", index_text), ("executive-summary.html", summary_text)):
            check(c["filename"] in text and "download" in text,
                  f"current download missing from {name}")
        check(c["tracking_name"] in index_text + summary_text, "tracking name missing from the site")
        check(c["note"] in summary_text, "submission note missing from the submission guide")
        root_page = (ROOT / "index.html").read_text()
        check("docs/index.html" in root_page and c["filename"] in root_page,
              "root route is stale or missing the download")
        readme = (ROOT / "README.md").read_text()
        embedded = readme.split("<!-- BEGIN USER BRIEF 20261006 -->", 1)[1] \
                         .split("<!-- END USER BRIEF 20261006 -->", 1)[0].strip()
        check(embedded == (DOCS / "user-prompt-20261006.md").read_text().strip(),
              "README does not retain the complete current brief")
        check("Maximize P(Win)" in readme and "Own the Outcome" in readme, "core values lost")
        check(c["filename"] in readme, "README does not point at the current download")

        # retained history: earlier artifacts keep their own exact bytes and no-go
        for name, digest_expected in {
            "docs/downloads/gemsdoe40-h4-contact-offset-depthkde-20261006-fab9f6619c02.tif":
                "ee73ffd76fabbaa1a2e77e17f57947a7db858916d713801e0c3502e49f6acabb",
            "docs/downloads/gemsdoe40-h7-rtp-euler-gravity-context-3d-kde-20261006-998f660f.tif":
                "e229aa9af26018bc80b0cca880d484a5b8362c26ad79ca362423953c32e66597",
            "docs/downloads/gemsdoe40-h2b-tmi-euler-natural-support-20261005.tif":
                "02486eaa491da2d8ebc6b1cd53ed2bdc8e995f316be1110a64eec484779b5a68",
        }.items():
            path = ROOT / name
            if path.is_file():
                check(digest(path) == digest_expected, f"retained artifact bytes changed: {name}")
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
    print(f"[site] OK: {len(list(DOCS.rglob('*.html')))} documentation pages plus root redirect; current H8 file/cloud bytes, novelty gate, registered settings, prior inventory, retained history and full prompt verified; no files rewritten")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
