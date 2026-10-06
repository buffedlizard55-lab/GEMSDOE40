#!/usr/bin/env python3
"""Non-destructive integrity check for the checked-in GitHub Pages site.

The pre-merge main-branch generator hard-coded the H40-4 artifact as "UPLOAD THIS" and wrote
several root pages and assets unconditionally. That advice predates this branch's expanded prior
corpus and H2-B no-go result. It is intentionally not allowed to overwrite the current site.
The site is now versioned as static HTML/JSON; this command validates local links and evidence
consistency without modifying files.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlsplit

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
CANDIDATE = DOCS / "downloads" / "gemsdoe40-h2b-tmi-euler-natural-support-20261005.tif"
VALIDATION = DOCS / "data" / "validation-h2b-20261005.json"
INVENTORY = DOCS / "data" / "prior_raster_inventory.json"
UNIQUENESS = DOCS / "data" / "uniqueness_audit.json"


class LocalLinks(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.targets: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        for key, value in attrs:
            if key in {"href", "src"} and value:
                self.targets.append(value)


def verify_pages() -> list[str]:
    problems: list[str] = []
    pages = [ROOT / "index.html", *sorted(DOCS.rglob("*.html"))]
    for page in pages:
        parser = LocalLinks()
        parser.feed(page.read_text(encoding="utf-8"))
        for target in parser.targets:
            parts = urlsplit(target)
            if parts.scheme or parts.netloc or not parts.path:
                continue
            local = (page.parent / parts.path).resolve()
            if ROOT not in local.parents and local != ROOT:
                problems.append(f"link escapes repository: {page.relative_to(ROOT)} -> {target}")
            elif not local.exists():
                problems.append(f"missing local link: {page.relative_to(ROOT)} -> {target}")
    return problems


def verify_evidence() -> list[str]:
    problems: list[str] = []
    root_page = (ROOT / "index.html").read_text(encoding="utf-8")
    if "HOLD — DO NOT SUBMIT" not in root_page or "docs/index.html" not in root_page:
        problems.append("root landing page must redirect to the current HOLD status page")
    if "UPLOAD THIS" in root_page or "Recommended submission" in root_page:
        problems.append("root landing page contains superseded upload advice")
    try:
        report = json.loads(VALIDATION.read_text(encoding="utf-8"))
        inventory = json.loads(INVENTORY.read_text(encoding="utf-8"))
        uniqueness = json.loads(UNIQUENESS.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return [f"cannot load current audit evidence: {exc}"]

    expected_sha = "02486eaa491da2d8ebc6b1cd53ed2bdc8e995f316be1110a64eec484779b5a68"
    actual_sha = hashlib.sha256(CANDIDATE.read_bytes()).hexdigest() if CANDIDATE.exists() else None
    if actual_sha != expected_sha or report.get("candidate_sha256") != expected_sha:
        problems.append(f"candidate SHA-256 mismatch: file={actual_sha}, report={report.get('candidate_sha256')}")
    gate = report.get("promotion_gate", {})
    if gate.get("final_action") != "HOLD; DO NOT SUBMIT" or gate.get("passed") is not False:
        problems.append("current candidate decision is not the expected explicit no-go")
    inv_count = len(inventory.get("unique_rasters", []))
    audit_count = report.get("prior_baselines", {}).get("inventory_unique_rasters")
    if inv_count != 279 or audit_count != inv_count:
        problems.append(f"prior inventory count mismatch: inventory={inv_count}, report={audit_count}")
    if report.get("prior_baselines", {}).get("same_grid_scored") != 277:
        problems.append("current proxy validation does not account for all 277 exact-grid priors")
    if uniqueness.get("same_grid_comparisons") != 277 or uniqueness.get("near_duplicate_count") != 0 or not uniqueness.get("uniqueness_pass"):
        problems.append("expanded H2-B uniqueness audit is missing, incomplete, or failing")
    if len(report.get("prior_raster_scores", [])) != inv_count:
        problems.append("validation report prior-score row count does not match inventory")
    return problems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="check static pages and evidence (default)")
    parser.parse_args()
    errors = verify_pages() + verify_evidence()
    if errors:
        print("[site] integrity check failed:")
        for error in errors:
            print(f"  - {error}")
        return 1
    print(f"[site] OK: {len(list(DOCS.rglob('*.html')))} documentation pages plus root redirect, "
          "277 uniqueness comparisons, 279 pinned prior blobs, and explicit HOLD status; no files rewritten")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
