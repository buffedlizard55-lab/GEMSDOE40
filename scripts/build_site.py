#!/usr/bin/env python3
"""Non-destructive local-link, release-byte and scientific-evidence checks.

Successor to the session-2 checker: the current release is the H8 lineament-weighted
cross-family SI = 0 depth-clustering candidate, whose chain is verified here, and every
earlier artifact (H13, H8 trace-locked, H8-ASA, H7, H4, H40, H2-B) keeps its own frozen
byte identity and its own explicit no-go status. Use ``scripts/build_h8_site.py`` (wrapped
by ``scripts/build_contact_site.py``) to rebuild the HTML from receipts. This checker never
rewrites the site and never turns a proxy failure into upload advice.
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
    """Prefer the recorded path; fall back to the published download copy."""
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
        session2 = load("session2-artifacts-20261006.json")
        audit, generation = load("h8-lineament-audit.json"), load("h8-lineament-generation.json")

        # ---------------------------------------------------------------- current release --
        check(c["status"] == "BUILT, AUDITED, NOT PROMOTED" and c["slot_eligible"] is False
              and c["weekly_submission_used"] is False and c["organizer_score"] is None,
              "current release status became an unverified claim")
        candidate = ROOT / c["path"]
        check(candidate.is_file() and digest(candidate) == c["sha256"], "current download SHA-256 mismatch")
        check(candidate.name == c["filename"] and candidate.stat().st_size == c["bytes"],
              "current file name/byte-size mismatch")
        check(c["canonical_pixels_sha256"] == audit["candidate"]["canonical_pixels_sha256"],
              "canonical pixel identity disagrees with the audit receipt")
        check(audit["published"]["sha256"] == c["sha256"], "publication and audit disagree on the file hash")
        check(audit["candidate"]["valid"] is True, "current file is not format-verified")
        check(audit["emitted_dots"] == c["dots"] == audit["candidate"]["candidate"]["in_footprint_nonzero"],
              "emitted dot count disagrees between receipt and audit")
        check(abs(c["mass"] - audit["mass"]) < 1e-6, "emitted mass disagrees between receipt and audit")
        check(audit["novelty"]["novel"] is True and not audit["novelty"]["exact_duplicates"],
              "novelty gate did not pass")
        check(audit["novelty"]["max_abs_pearson"] <= audit["novelty"]["thresholds"]["max_abs_pearson"],
              "raw correlation exceeds the registered novelty threshold")
        check(audit["novelty"]["max_jaccard_topmass"] <= audit["novelty"]["thresholds"]["max_jaccard_topmass"],
              "top-mass Jaccard exceeds the registered novelty threshold")
        check(generation["construction_reads_labels_or_priors"] is False,
              "generation claims to read labels/priors")
        check(generation["configuration"]["si"] == 0, "registered structural index changed after generation")
        check(abs(generation["configuration"]["nms_spacing_px"] - 2.8) < 1e-9,
              "registered emission spacing changed after generation")
        check(digest(DOCS / "downloads" / "h8-lineament-solutions.csv.gz") == c["cloud"]["sha256"],
              "downloadable depth-labelled solution cloud changed")
        check(c["note_characters"] == len(c["note"]) <= 200, "submission note length is wrong")
        check(digest(ROOT / c["hard_twin"]["path"]) == c["hard_twin"]["sha256"], "hard twin changed")

        # ------------------------------------------- retained session-2/3 chains (unchanged) --
        h13g, h13f, h13u, h13v = (load(n) for n in (
            "h13-generation.json", "h13-format.json", "h13-uniqueness.json", "h13-validation.json"))
        h8g, h8f, h8u, h8v = (load(n) for n in (
            "h8-generation.json", "h8-format.json", "h8-uniqueness.json", "h8-validation.json"))
        h4g, h4u, h4v, h4repro = (load(n) for n in (
            "h4-generation.json", "h4-uniqueness.json", "h4-validation.json", "h4-reproduction.json"))

        def verify_chain(tag, gen, fmt, uniq, val, path_key="path"):
            artifact = resolve_artifact(gen["candidate_format"][path_key])
            check(artifact.is_file() and digest(artifact) == gen["candidate_format"]["sha256"],
                  f"{tag}: download SHA-256 mismatch")
            check(gen["candidate_format"]["sha256"] == fmt["sha256"] == uniq["candidate_sha256"]
                  == val["candidate"]["sha256"], f"{tag}: evidence disagrees on candidate identity")
            check(fmt["valid"] is True and gen["candidate_format"]["valid"] is True,
                  f"{tag}: file is not format-verified")
            check(fmt["canonical_pixels_sha256"] == uniq["candidate_canonical_sha256"]
                  == val["candidate"]["canonical_pixels_sha256"], f"{tag}: canonical pixel identities disagree")
            n = len(inv["unique_rasters"])
            check(n == inv["unique_git_blobs"] == uniq["inventory_blobs"] == uniq["hashed_and_audited_blobs"]
                  == len(uniq["comparisons"]), f"{tag}: incomplete raw-output audit counts")
            check(uniq["inventory_sha256"] == digest(DOCS / "data/prior-inventory-20261006.json"),
                  f"{tag}: prior inventory changed since the novelty audit")
            check(uniq["completeness_pass"] and uniq["uniqueness_pass"] and uniq["near_duplicate_count"] == 0
                  and not uniq["issues"], f"{tag}: novelty/integrity does not pass")
            check(val["gate"]["all_pass"] is False and val["gate"]["g3_pass"] is False,
                  f"{tag}: frozen proxy gate unexpectedly changed to pass")
            check(val["gate"]["weekly_slot_used"] is False and val["gate"]["organizer_score"] is None,
                  f"{tag}: unverified submission/score claim")
            check(gen["construction_reads_proxy_or_prior_predictions"] is False,
                  f"{tag}: generation depends on old predictions/proxy")
            prereg_doc = ROOT / gen["preregistration"]["path"]
            check(prereg_doc.is_file(), f"{tag}: preregistration missing")
            check(val["preregistration"]["sha256_at_generation"] == gen["preregistration"]["sha256"],
                  f"{tag}: generation-time preregistration hash not reproduced")
            if val["preregistration"]["amended_after_generation"]:
                check("amended after generation" in val["preregistration"]["note"],
                      f"{tag}: post-generation amendment not disclosed")
            else:
                check(digest(prereg_doc) == gen["preregistration"]["sha256"],
                      f"{tag}: frozen registration changed after generation")
            for source, expected in gen["code_sha256"].items():
                check(digest(ROOT / source) == expected, f"{tag}: scientific source changed after run: {source}")

        check(session2["experiment"] == "H13" and session2["slot_eligible"] is False
              and session2["status"] == "HOLD — DO NOT SUBMIT", "retained session-2 no-go changed")
        check(session2["organizer_score"] is None, "unverified organizer-score claim in retained manifest")
        retained = ROOT / session2["path"]
        check(retained.is_file() and digest(retained) == session2["sha256"],
              "retained H13 download SHA-256 mismatch")
        check(session2["canonical_pixels_sha256"] == h13f["canonical_pixels_sha256"],
              "retained canonical pixel identity mismatch")
        check(session2["note_characters"] == len(session2["note"]) <= 200, "retained note length wrong")
        check(session2["cloud"]["sha256"] == digest(ROOT / session2["cloud"]["path"]),
              "retained H4 solution cloud has changed")
        verify_chain("H13", h13g, h13f, h13u, h13v)
        verify_chain("H8-tracelock", h8g, h8f, h8u, h8v)

        h8asag, h8asaf, h8asau, h8asav, h8asarepro = (load(n) for n in (
            "h8asa-generation.json", "h8asa-format.json", "h8asa-uniqueness.json",
            "h8asa-validation.json", "h8asa-reproduction.json"))
        hc = session2["h8asa_candidate"]
        h8asa_path = ROOT / hc["path"]
        check(h8asa_path.is_file() and digest(h8asa_path) == hc["sha256"], "H8-ASA download SHA-256 mismatch")
        check(hc["sha256"] == h8asaf["sha256"] == h8asag["candidate_format"]["sha256"]
              == h8asau["candidate_sha256"] == h8asav["candidate_sha256"],
              "H8-ASA evidence disagrees on candidate identity")
        check(h8asaf["valid"] is True and h8asag["candidate_format"]["valid"] is True,
              "H8-ASA file is not format-verified")
        check(hc["canonical_pixels_sha256"] == h8asaf["canonical_pixels_sha256"]
              == h8asau["candidate_canonical_sha256"], "H8-ASA canonical pixel identities disagree")
        check(len(inv["unique_rasters"]) == h8asau["inventory_blobs"] == h8asau["hashed_and_audited_blobs"]
              == len(h8asau["comparisons"]), "H8-ASA incomplete audit counts")
        check(h8asau["completeness_pass"] and h8asau["uniqueness_pass"] and h8asau["near_duplicate_count"] == 0
              and not h8asau["issues"], "H8-ASA novelty/integrity does not pass")
        check(h8asav["promotion_gate"]["passed"] is False and h8asav["promotion_gate"]["slot_eligible"] is False
              and h8asav["promotion_gate"]["organizer_score"] is None
              and h8asav["promotion_gate"]["weekly_submission_used"] is False, "H8-ASA no-go changed")
        check(hc["status"] == h8asav["promotion_gate"]["final_action"] == "HOLD — DO NOT SUBMIT",
              "H8-ASA no-go is not explicit")
        check(hc["portal_safe_twin"]["sha256"] == digest(DOCS / "downloads" / hc["portal_safe_twin"]["filename"]),
              "H8-ASA portal twin changed")
        check(hc["zip"]["sha256"] == digest(DOCS / "downloads" / hc["zip"]["filename"]), "H8-ASA zip changed")
        check(hc["cloud"]["sha256"] == digest(ROOT / hc["cloud"]["path"]), "H8-ASA cloud changed")
        check(hc["note_characters"] == len(hc["note"]) <= 200, "H8-ASA note length wrong")
        check(h8asarepro["pass"] and h8asarepro["actual_tiff_sha256"] == hc["sha256"]
              and h8asarepro["actual_cloud_sha256"] == hc["cloud"]["sha256"]
              and h8asarepro["actual_zeros_sha256"] == hc["portal_safe_twin"]["sha256"]
              and h8asarepro["actual_zip_sha256"] == hc["zip"]["sha256"], "H8-ASA byte reproduction disagrees")
        check(h8asag["construction_uses_proxy_or_prior_predictions"] is False,
              "H8-ASA generation depends on old predictions/proxy")
        check(digest(DOCS / "research/h8asa-preregistration-20261006.md") == h8asag["preregistration_sha256"],
              "H8-ASA frozen registration changed after generation")
        for source, expected in h8asag["code_sha256"].items():
            check(digest(ROOT / source) == expected, f"H8-ASA scientific source changed after run: {source}")
        sib = session2["sibling_candidate"]
        check(digest(ROOT / sib["path"]) == sib["sha256"] and sib["status"].startswith("HOLD"),
              "H8 sibling record mismatch")
        arch = session2["archived_candidate_h4"]
        check(digest(ROOT / "docs/downloads" / arch["filename"]) == arch["sha256"] and arch["status"].startswith("HOLD"),
              "archived H4 record mismatch")
        check(h4v["promotion_gate"]["final_action"] == "HOLD — DO NOT SUBMIT"
              and h4v["promotion_gate"]["passed"] is False, "archived H4 no-go changed")
        check(h4repro["pass"] and h4repro["actual_tiff_sha256"] == arch["sha256"]
              and h4repro["actual_cloud_sha256"] == session2["cloud"]["sha256"],
              "archived H4 byte reproduction disagrees")

        # ------------------------------------------------------------------ shared records --
        projects_check = (projects["project_count"] == len(projects["projects"]) == 44)
        check(projects_check, "complete 44-project source register missing")
        check(projects["reported_score_count"]
              == sum(s["score"] is not None for p in projects["projects"] for s in p["submissions"]) == 48,
              "48 owner-reported scores were not preserved")
        cached = len(list((ROOT / "data" / "prior").glob("*.tif")))
        comparable = audit["novelty"]["corpus_size"] + audit["novelty"]["unreadable_or_incomparable"]
        check(inv["unique_git_blobs"] == len(inv["unique_rasters"]) >= cached,
              "prior inventory is smaller than the cached corpus")
        check(cached == comparable, "cached prior corpus disagrees with the novelty audit")

        index_text = (DOCS / "index.html").read_text()
        summary_text = (DOCS / "executive-summary.html").read_text()
        for name, text in (("index.html", index_text), ("executive-summary.html", summary_text)):
            check(c["filename"] in text and "download" in text and "HOLD — DO NOT SUBMIT" in text,
                  f"current download/no-go missing from {name}")
            check(session2["filename"] in text and sib["filename"] in text and hc["filename"] in text,
                  f"retained session-2 downloads missing from {name}")
        check("HOLD — DO NOT SUBMIT H40" in index_text and "HOLD — DO NOT SUBMIT H7 or H40" in index_text,
              "retained HOLD statements missing from the index")
        check(c["tracking_name"] in index_text + summary_text, "tracking name missing from the site")
        check(c["note"] in summary_text, "submission note missing from the submission guide")
        root_page = (ROOT / "index.html").read_text()
        check("docs/index.html" in root_page and c["filename"] in root_page
              and "HOLD — DO NOT SUBMIT" in root_page, "root route is stale or missing the download")
        check("UPLOAD THIS" not in root_page and "Recommended submission" not in root_page,
              "root still promotes archived output")
        readme = (ROOT / "README.md").read_text()
        embedded = readme.split("<!-- BEGIN USER BRIEF 20261006 -->", 1)[1] \
                         .split("<!-- END USER BRIEF 20261006 -->", 1)[0].strip()
        check(embedded == (DOCS / "user-prompt-20261006.md").read_text().strip(),
              "README does not retain the complete current brief")
        check("Maximize P(Win)" in readme and "Own the Outcome" in readme, "core values lost")
        check(c["filename"] in readme, "README does not point at the current download")

        h7 = load("validation-h7-20261006.json")
        h7_path = DOCS / "downloads/gemsdoe40-h7-rtp-euler-gravity-context-3d-kde-20261006-998f660f.tif"
        check(digest(h7_path) == h7["candidate_format_receipt"]["sha256"], "retained H7 bytes changed")
        check(h7["slot_eligible"] is False and h7["weekly_slot_used"] is False and h7["organizer_score"] is None,
              "retained H7 no-go changed")
        h40 = json.loads((DOCS / "downloads/gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-audit.json").read_text())
        for suffix, key in (("zeros", "zeros_tif"), ("nan", "nan_tif")):
            path = DOCS / f"downloads/gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-{suffix}.tif"
            check(digest(path) == h40[key]["sha256"], f"retained H40 {suffix} bytes changed")
        check((DOCS / "reports/integration-20261006.md").is_file(),
              "parallel experiment/proxy-profile integration notes missing")
        old = DOCS / "downloads/gemsdoe40-h2b-tmi-euler-natural-support-20261005.tif"
        check(digest(old) == "02486eaa491da2d8ebc6b1cd53ed2bdc8e995f316be1110a64eec484779b5a68",
              "historical H2-B bytes changed")
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
    print(f"[site] OK: {len(list(DOCS.rglob('*.html')))} documentation pages plus root redirect; current H8 "
          "lineament file/cloud/twin bytes, novelty gate, registered settings, retained H13/H8/H8-ASA/H4/H7/H40 "
          "identities, complete prior inventory, frozen no-go and full prompt verified; no files rewritten")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
