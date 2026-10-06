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


def resolve_artifact(path_field: str) -> Path:
    """Prefer the recorded path; fall back to the published download copy.

    Generation receipts record the transient ``work/`` artifact (gitignored by
    design). The byte-identical published download lives in ``docs/downloads/``;
    the sha chain across generation/format/uniqueness/validation receipts proves
    identity, so either copy validates the same digest.
    """
    primary = ROOT / path_field
    if primary.is_file():
        return primary
    fallback = DOCS / "downloads" / Path(path_field).name
    return fallback if fallback.is_file() else primary


def verify_evidence() -> list[str]:
    problems = []
    def check(condition, message):
        if not condition: problems.append(message)
    def load(name): return json.loads((DOCS / "data" / name).read_text())
    try:
        c, inv, projects = (load(n) for n in (
            "current-candidate.json", "prior-inventory-20261006.json", "project-site-review-20261006.json"))
        h13g, h13f, h13u, h13v = (load(n) for n in (
            "h13-generation.json", "h13-format.json", "h13-uniqueness.json", "h13-validation.json"))
        h8g, h8f, h8u, h8v = (load(n) for n in (
            "h8-generation.json", "h8-format.json", "h8-uniqueness.json", "h8-validation.json"))
        h4g, h4u, h4v, h4repro = (load(n) for n in (
            "h4-generation.json", "h4-uniqueness.json", "h4-validation.json", "h4-reproduction.json"))

        def verify_chain(tag, gen, fmt, uniq, val, path_key="path"):
            candidate = resolve_artifact(gen["candidate_format"][path_key])
            check(candidate.is_file() and digest(candidate) == gen["candidate_format"]["sha256"], f"{tag}: download SHA-256 mismatch")
            check(gen["candidate_format"]["sha256"] == fmt["sha256"] == uniq["candidate_sha256"] == val["candidate"]["sha256"], f"{tag}: evidence disagrees on candidate identity")
            check(fmt["valid"] is True and gen["candidate_format"]["valid"] is True, f"{tag}: file is not format-verified")
            check(fmt["canonical_pixels_sha256"] == uniq["candidate_canonical_sha256"] == val["candidate"]["canonical_pixels_sha256"], f"{tag}: canonical pixel identities disagree")
            n = len(inv["unique_rasters"])
            check(n == inv["unique_git_blobs"] == uniq["inventory_blobs"] == uniq["hashed_and_audited_blobs"] == len(uniq["comparisons"]), f"{tag}: incomplete raw-output audit counts")
            check(uniq["inventory_sha256"] == digest(DOCS / "data/prior-inventory-20261006.json"), f"{tag}: prior inventory changed since the novelty audit")
            check(uniq["completeness_pass"] and uniq["uniqueness_pass"] and uniq["near_duplicate_count"] == 0 and not uniq["issues"], f"{tag}: novelty/integrity does not pass")
            check(val["gate"]["all_pass"] is False and val["gate"]["g3_pass"] is False, f"{tag}: frozen proxy gate unexpectedly changed to pass")
            check(val["gate"]["weekly_slot_used"] is False and val["gate"]["organizer_score"] is None, f"{tag}: unverified submission/score claim")
            check(gen["construction_reads_proxy_or_prior_predictions"] is False, f"{tag}: generation depends on old predictions/proxy")
            prereg_doc = ROOT / gen["preregistration"]["path"]
            check(prereg_doc.is_file(), f"{tag}: preregistration missing")
            check(val["preregistration"]["sha256_at_generation"] == gen["preregistration"]["sha256"], f"{tag}: generation-time preregistration hash not reproduced")
            if val["preregistration"]["amended_after_generation"]:
                check("amended after generation" in val["preregistration"]["note"], f"{tag}: post-generation amendment not disclosed")
            else:
                check(digest(prereg_doc) == gen["preregistration"]["sha256"], f"{tag}: frozen registration changed after generation")
            for source, expected in gen["code_sha256"].items():
                check(digest(ROOT / source) == expected, f"{tag}: scientific source changed after run: {source}")

        # Current candidate: H13 (session 2), sibling H8, archived H4.
        check(c["experiment"] == "H13" and c["slot_eligible"] is False and c["status"] == "HOLD — DO NOT SUBMIT", "current candidate no-go changed")
        check(c["organizer_score"] is None, "unverified organizer-score claim")
        candidate = ROOT / c["path"]
        check(candidate.is_file() and digest(candidate) == c["sha256"], "current download SHA-256 mismatch")
        check(c["filename"] == candidate.name and candidate.stat().st_size == c["bytes"], "current file name/byte-size mismatch")
        check(c["canonical_pixels_sha256"] == h13f["canonical_pixels_sha256"], "current canonical pixel identity mismatch")
        check(c["note_characters"] == len(c["note"]) <= 200, "submission note exceeds portal limit or has wrong recorded length")
        check(c["cloud"]["sha256"] == digest(ROOT / c["cloud"]["path"]), "downloadable solution cloud has changed")
        verify_chain("H13", h13g, h13f, h13u, h13v)
        verify_chain("H8", h8g, h8f, h8u, h8v)
        sib = c["sibling_candidate"]
        check(digest(ROOT / sib["path"]) == sib["sha256"] and sib["status"].startswith("HOLD"), "H8 sibling record mismatch")
        arch = c["archived_candidate_h4"]
        check(digest(ROOT / "docs/downloads" / arch["filename"]) == arch["sha256"] and arch["status"].startswith("HOLD"), "archived H4 record mismatch")
        check(h4v["promotion_gate"]["final_action"] == "HOLD — DO NOT SUBMIT" and h4v["promotion_gate"]["passed"] is False, "archived H4 no-go changed")
        check(h4repro["pass"] and h4repro["actual_tiff_sha256"] == arch["sha256"] and h4repro["actual_cloud_sha256"] == c["cloud"]["sha256"], "archived H4 byte reproduction disagrees")
        check(projects["project_count"] == len(projects["projects"]) == 44, "complete 44-project source register missing")
        check(projects["reported_score_count"] == sum(s["score"] is not None for p in projects["projects"] for s in p["submissions"]) == 48, "48 owner-reported scores were not preserved")
        for name in ("index.html", "executive-summary.html"):
            text = (DOCS / name).read_text()
            check(c["filename"] in text and "download" in text.lower() and "HOLD — DO NOT SUBMIT" in text, f"current download/no-go missing from {name}")
            check(c["sibling_candidate"]["filename"] in text, f"H8 sibling download missing from {name}")
        root_page = (ROOT / "index.html").read_text()
        check("docs/index.html" in root_page and c["filename"] in root_page and "HOLD — DO NOT SUBMIT" in root_page, "root route is stale or missing the download")
        check("UPLOAD THIS" not in root_page and "Recommended submission" not in root_page, "root still promotes archived output")
        readme = (ROOT / "README.md").read_text()
        embedded = readme.split("<!-- BEGIN USER BRIEF 20261006 -->", 1)[1].split("<!-- END USER BRIEF 20261006 -->", 1)[0].strip()
        check(embedded == (DOCS / "user-prompt-20261006.md").read_text().strip(), "README does not retain the complete current brief")
        check("Maximize P(Win)" in readme and "Own the Outcome" in readme, "core values lost")
        # Concurrent mainline artifacts retain their own exact identities and no-go.
        h7 = load("validation-h7-20261006.json")
        h7_path = DOCS / "downloads/gemsdoe40-h7-rtp-euler-gravity-context-3d-kde-20261006-998f660f.tif"
        check(digest(h7_path) == h7["candidate_format_receipt"]["sha256"], "retained H7 bytes changed")
        check(h7["slot_eligible"] is False and h7["weekly_slot_used"] is False and h7["organizer_score"] is None, "retained H7 no-go changed")
        h40 = json.loads((DOCS / "downloads/gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-audit.json").read_text())
        for suffix, key in (("zeros", "zeros_tif"), ("nan", "nan_tif")):
            path = DOCS / f"downloads/gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-{suffix}.tif"
            check(digest(path) == h40[key]["sha256"], f"retained H40 {suffix} bytes changed")
        check((DOCS / "reports/integration-20261006.md").is_file(), "parallel experiment/proxy-profile integration notes missing")
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
