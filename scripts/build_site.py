#!/usr/bin/env python3
"""Non-destructive integrity check for the checked-in GitHub Pages site.

The original main-branch generator hard-coded an obsolete upload recommendation. This script
validates the current H7 decision and local links without rewriting versioned HTML/JSON.
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
CANDIDATE = DOCS / "downloads" / "gemsdoe40-h7-rtp-euler-gravity-context-3d-kde-20261006-998f660f.tif"
VALIDATION = DOCS / "data" / "validation-h7-20261006.json"
INVENTORY = DOCS / "data" / "prior_raster_inventory-20261006.json"
UNIQUENESS = DOCS / "data" / "uniqueness-audit-h7-20261006.json"
PRIOR_SCORES = DOCS / "data" / "prior-holdout-scores-h7-20261006.json"
#: second current artifact (2026-10-06, H40 Euler depth-cluster) and its receipt
CANDIDATE_H40 = DOCS / "downloads" / "gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-zeros.tif"
RECEIPT_H40 = DOCS / "downloads" / "gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-audit.json"
SHA_H40 = "57896abee36d6722f587f543e1f65ddfd9af19163e20653e852bfb3162b10902"


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
    site_home = (DOCS / "index.html").read_text(encoding="utf-8")
    executive = (DOCS / "executive-summary.html").read_text(encoding="utf-8")
    if "HOLD — DO NOT SUBMIT" not in root_page or "docs/index.html" not in root_page:
        problems.append("root landing page must redirect to the current HOLD status page")
    if "UPLOAD THIS" in root_page or "Recommended submission" in root_page:
        problems.append("root landing page contains superseded upload advice")
    if "H7 research-only GeoTIFF" not in site_home or "HOLD — DO NOT SUBMIT" not in site_home:
        problems.append("current site home must make the H7 research download and no-go status obvious")
    if "No candidate is cleared today" not in executive or "Do not upload H7" not in executive:
        problems.append("executive summary must retain the no-upload guide for the current research-only candidate")

    try:
        report = json.loads(VALIDATION.read_text(encoding="utf-8"))
        inventory = json.loads(INVENTORY.read_text(encoding="utf-8"))
        uniqueness = json.loads(UNIQUENESS.read_text(encoding="utf-8"))
        prior_scores = json.loads(PRIOR_SCORES.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return problems + [f"cannot load current audit evidence: {exc}"]

    expected_sha = "e229aa9af26018bc80b0cca880d484a5b8362c26ad79ca362423953c32e66597"
    actual_sha = hashlib.sha256(CANDIDATE.read_bytes()).hexdigest() if CANDIDATE.exists() else None
    receipt = report.get("candidate_format_receipt", {})
    if actual_sha != expected_sha or receipt.get("sha256") != expected_sha:
        problems.append(f"H7 candidate SHA-256 mismatch: file={actual_sha}, report={receipt.get('sha256')}")
    if receipt.get("valid") is not True:
        problems.append("H7 candidate format receipt is not a pass")
    gate = report.get("promotion_gate_vs_highest_all_prior", {})
    if gate.get("passed") is not False or report.get("slot_eligible") is not False:
        problems.append("H7 report does not preserve the expected explicit no-go")
    if report.get("weekly_slot_used") is not False or report.get("organizer_score") is not None:
        problems.append("H7 report must state that no slot was used and no organizer score exists")

    inventory_rows = inventory.get("unique_rasters", [])
    refresh = inventory.get("current_head_refresh", {})
    if len(inventory_rows) != 306 or inventory.get("unique_pinned_repo_ref_path_locations") != 541:
        problems.append("refreshed prior inventory count/location mismatch")
    if refresh.get("owner_public_repositories_matching_gemsdoe") != 55 or refresh.get("current_head_tiff_paths_scanned") != 501:
        problems.append("refreshed current-head repository/path count mismatch")
    if refresh.get("recursive_trees_truncated") != 0 or refresh.get("unclassified_tiff_paths_excluded_pending_manual_review"):
        problems.append("prior inventory is truncated or has pending row-level classification")

    prior_summary = report.get("prior_corpus_holdout", {})
    if prior_summary.get("same_grid_prior_outputs_scored") != 304 or prior_summary.get("prior_format_eligible_count") != 183:
        problems.append("current-proxy prior comparison does not account for the 304 same-grid outputs and 183 locally format-eligible subset")
    if len(prior_scores.get("scores", [])) != 304 or prior_scores.get("proxy_sha256") != report.get("input_sha256", {}).get("current_proxy"):
        problems.append("complete same-grid prior score receipt is missing or uses a different proxy pin")
    if uniqueness.get("inventory_unique_blobs") != 306 or uniqueness.get("same_grid_comparisons") != 304:
        problems.append("H7 uniqueness audit does not cover all 304 exact-grid prior blobs")
    if uniqueness.get("missing_prior_cache_count") != 0 or uniqueness.get("near_duplicate_count") != 0 or not uniqueness.get("uniqueness_pass"):
        problems.append("H7 uniqueness audit is incomplete or failing")
    if uniqueness.get("top_budget_for_comparison") != 45_962:
        problems.append("H7 uniqueness audit did not use the frozen 45,962-cell top budget")
    # ---- second current artifact: H40 Euler depth-cluster ------------------------------------
    if not CANDIDATE_H40.exists():
        problems.append(f"H40 artifact missing: {CANDIDATE_H40.name}")
    else:
        got40 = hashlib.sha256(CANDIDATE_H40.read_bytes()).hexdigest()
        if got40 != SHA_H40:
            problems.append(f"H40 artifact SHA-256 mismatch: {got40}")
    try:
        receipt40 = json.loads(RECEIPT_H40.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        problems.append(f"H40 receipt unreadable: {exc}")
    else:
        if receipt40.get("zeros_tif", {}).get("sha256") != SHA_H40:
            problems.append("H40 receipt does not pin the artifact hash")
        inst = receipt40.get("stage", {}).get("emission", {}).get("instrument", {})
        pred = inst.get("predicted_live")
        if pred is None or pred >= 0.2778:
            problems.append(f"H40 artifact is not held: instrument prediction {pred}")
        if inst.get("lm_calibrated") is None or inst["lm_calibrated"] >= (inst.get("lm_incumbent") or 0.0):
            problems.append("H40 LM instrument does not show a shortfall against its incumbent")
        if receipt40.get("stage", {}).get("uniqueness", {}).get("is_new") is not True:
            problems.append("H40 uniqueness audit did not pass")
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
    print(f"[site] OK: {len(list(DOCS.rglob('*.html')))} documentation pages plus root redirect; "
          "H7 HOLD decision, complete 306-blob inventory, 304 holdout/uniqueness comparisons, and no upload advice")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
