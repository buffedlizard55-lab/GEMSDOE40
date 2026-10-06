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
        c, g, f, u, v, inv, repro, projects = (load(n) for n in (
            "current-candidate.json", "h4-generation.json", "h4-format.json", "h4-uniqueness.json",
            "h4-validation.json", "prior-inventory-20261006.json", "h4-reproduction.json", "project-site-review-20261006.json"))
        candidate = ROOT / c["path"]
        check(candidate.is_file() and digest(candidate) == c["sha256"], "current download SHA-256 mismatch")
        check(c["filename"] == candidate.name and candidate.stat().st_size == c["bytes"], "current file name/byte-size mismatch")
        check(c["sha256"] == f["sha256"] == g["candidate_format"]["sha256"] == u["candidate_sha256"] == v["candidate_sha256"], "current evidence disagrees on candidate identity")
        check(f["valid"] is True and g["candidate_format"]["valid"] is True, "current file is not format-verified")
        check(c["cloud"]["sha256"] == digest(ROOT / c["cloud"]["path"]), "downloadable solution cloud has changed")
        check(c["note_characters"] == len(c["note"]) <= 200, "submission note exceeds portal limit or has wrong recorded length")
        check(c["canonical_pixels_sha256"] == f["canonical_pixels_sha256"] == u["candidate_canonical_sha256"], "canonical pixel identities disagree")
        n = len(inv["unique_rasters"])
        check(n == inv["unique_git_blobs"] == u["inventory_blobs"] == u["hashed_and_audited_blobs"] == len(u["comparisons"]), "incomplete raw-output audit counts")
        check(u["inventory_sha256"] == digest(DOCS / "data/prior-inventory-20261006.json"), "prior inventory changed since the novelty audit")
        check(u["completeness_pass"] and u["uniqueness_pass"] and u["near_duplicate_count"] == 0 and not u["issues"], "novelty/integrity does not pass")
        gate = v["promotion_gate"]
        check(c["slot_eligible"] is False and gate["slot_eligible"] is False and gate["passed"] is False, "current no-go unexpectedly changed; manually review publication advice")
        check(c["status"] == gate["final_action"] == "HOLD — DO NOT SUBMIT", "current no-go is not explicit")
        check(c["organizer_score"] is None and gate["organizer_score"] is None and gate["weekly_submission_used"] is False, "unverified organizer-score or submission claim")
        check(gate["blocks_with_truth"] == 16 and gate["empty_blocks"] == 8 and gate["blocks_required"] == 18 and not gate["strict_win_requirement_feasible"], "infeasible inherited gate was hidden or silently relaxed")
        check(g["construction_uses_proxy_or_prior_predictions"] is False, "current generation depends on old predictions/proxy")
        check(repro["pass"] and repro["actual_tiff_sha256"] == c["sha256"] and repro["actual_cloud_sha256"] == c["cloud"]["sha256"], "independent byte reproduction is absent or disagrees")
        prereg = digest(DOCS / "research/h4-preregistration-20261006.md")
        check(g["preregistration_sha256"] == prereg, "frozen registration changed after generation")
        for source, expected in {**g["code_sha256"], **v.get("audit_code_sha256", {})}.items():
            check(digest(ROOT / source) == expected, f"scientific source changed after run/audit: {source}")
        check(projects["project_count"] == len(projects["projects"]) == 44, "complete 44-project source register missing")
        check(projects["reported_score_count"] == sum(s["score"] is not None for p in projects["projects"] for s in p["submissions"]) == 48, "48 owner-reported scores were not preserved")
        for name in ("index.html", "executive-summary.html"):
            text = (DOCS / name).read_text()
            check(c["filename"] in text and "download" in text and "HOLD — DO NOT SUBMIT" in text, f"current download/no-go missing from {name}")
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
