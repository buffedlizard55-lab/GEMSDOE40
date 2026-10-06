#!/usr/bin/env python3
"""Refresh the prior-prediction TIFF corpus from all public owner repository heads.

Every owner/repository/default-branch head and recursive Git tree is recorded. TIFFs are explicitly
classified as prior output, proxy-circular output, input/template/truth, fixture/source, or review
needed. Only prediction/output rasters enter the uniqueness and holdout corpus; source maps,
proxy truth, labels, templates, and test fixtures do not. New output blobs are fetched through `gh`
and byte-verified against Git blob SHA-1 before their SHA-256 is admitted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
from urllib.parse import quote

OWNER = "buffedlizard55-lab"


def gh_json(endpoint: str) -> dict | list:
    return json.loads(subprocess.check_output(["gh", "api", endpoint], text=True))


def gh_raw(endpoint: str) -> bytes:
    return subprocess.check_output(["gh", "api", "-H", "Accept: application/vnd.github.raw", endpoint])


def git_blob_sha1(data: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def classify_tiff(path: str) -> tuple[str, bool]:
    """Return (evidence class, include in prior-output corpus)."""
    p = path.lower()
    base = p.rsplit("/", 1)[-1]

    # Explicitly proxy-trained/SGMC-derived output maps are retained for novelty checks and
    # separately flagged as circular if used against the matching proxy.
    if "proxy_only" in p or "sgmc-gap" in p or "sgmc_gap" in p:
        return "proxy-circular research/prediction output; include, flag as circular", True

    # Known raster inputs and source/context rasters. These can be scientifically relevant inputs,
    # but are not prior predictions and must not enter the prediction-to-prediction novelty corpus.
    if base in {
        "training_features.tif", "example_submission.tif", "sample_submission.tif",
        "existing_faults.tif", "labels.tif", "template-mask.tif", "fixture_labels.tif",
        "fixture_features_int16.tif", "derived_sgmc_faults_100m_u8.tif",
        "derived_gdr_2m_probes_100m_u8.tif", "derived_gdr_paleo_100m_u8.tif",
        "derived_gdr_qfaults_v2_100m_u8.tif", "derived_gdr_volcanics_100m_u8.tif",
        "proxy_catalogue.tif", "proxy_catalogue_sgmc.tif", "proxy_window.tif",
        "qfaults_catalogue.tif", "qfaults_prior_u8.tif", "lidar_scarp_features_u8.tif",
        "geodawn_extensions_u8.tif", "geodawn_rad_u8.tif", "eval_labels.tif", "train_labels.tif",
    }:
        return "input/template/label/proxy/source; exclude", False
    if any(token in p for token in (
        "/data/bridge/", "/legacy/data/bridge/", "/tests/fixture/", "/data/fixture/",
        "/external/audit_sources/", "/external/dem/", "/external/qfaults/",
        "/data/evidence/proxy/", "/data/evidence/xcat/", "/data/evidence/seghold/",
        "/external/geodawn_extensions/", "/external/geodawn_rad/",
    )):
        return "input/template/label/proxy/source/fixture; exclude", False

    # Explicit prior submissions, scored outputs, model predictions, research candidates and
    # archived/quarantined predictions are all needed for a conservative novelty audit.
    output_prefixes = (
        "docs/downloads/", "downloads/", "submissions/", "assets/lb_anchors/",
        "data/evidence/leaderboard_anchor/", "external/scored/", "inputs/",
        "archive/legacy_candidates/", "docs/archive/", "docs/research/quarantine/",
        "data/evidence/union_po_loo/", "data/evidence/newfault/", "data/evidence/suture/",
        "evidence/experiments/", "evidence/history/pre_correction_downloads/",
    )
    if p.startswith(output_prefixes) or "/downloads/" in p or "/submissions/" in p:
        return "prior prediction/submission/research output; include", True
    if p.startswith("docs/gemsdoe") and p.endswith((".tif", ".tiff")):
        return "prior prediction output at docs root; include", True
    if p.startswith("data/derived/context_detector_prob_"):
        return "model probability output (multiple variants); include", True
    if base in {"submission.tif", "submission_conformant.tif"}:
        return "prior submission/prediction output; include", True
    if any(token in base for token in ("prediction", "candidate", "prob_raw", "anchor")):
        return "likely prediction output; include after path/name context review", True

    return "unclassified TIFF; manual review required; exclude until classified", False


def _artifact_location(repo: str, ref: str, path: str, blob_sha: str, size: int) -> dict:
    return {
        "repo": repo,
        "ref": ref,
        "path": path,
        "bytes": int(size),
        "raw_url": f"https://raw.githubusercontent.com/{OWNER}/{repo}/{ref}/{quote(path, safe='/')}",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, default=Path("docs/data/prior_raster_inventory.json"))
    parser.add_argument("--out", type=Path, default=Path("docs/data/prior_raster_inventory-20261006.json"))
    parser.add_argument("--cache", type=Path, default=Path("/tmp/gemsdoe40-prior-cache"))
    args = parser.parse_args()
    base = json.loads(args.base.read_text(encoding="utf-8"))
    inventory = json.loads(json.dumps(base))
    entries = inventory["unique_rasters"]
    by_blob = {entry["git_blob_sha"]: entry for entry in entries}
    args.cache.mkdir(parents=True, exist_ok=True)

    repos = gh_json(f"users/{OWNER}/repos?per_page=100&type=public")
    repos = sorted((repo for repo in repos if "gemsdoe" in repo["name"].lower()), key=lambda r: r["name"])
    if len(repos) < 50:
        raise RuntimeError(f"owner repository listing unexpectedly short: {len(repos)} matches")

    heads: list[dict] = []
    scan_paths: list[dict] = []
    new_output_locations = 0
    new_unique_blobs = 0
    prediction_repositories: set[str] = set(inventory.get("source_repositories_with_tiff_outputs", []))
    excluded_review: list[dict] = []
    truncated: list[str] = []
    tree_count = 0

    for repository in repos:
        name = repository["name"]
        branch = repository["default_branch"]
        commit = gh_json(f"repos/{OWNER}/{name}/commits/{branch}")
        head = commit["sha"]
        tree = gh_json(f"repos/{OWNER}/{name}/git/trees/{head}?recursive=1")
        if tree.get("truncated"):
            truncated.append(name)
            continue
        tree_count += 1
        heads.append({"repo": name, "default_branch": branch, "head_commit": head})
        for item in tree.get("tree", []):
            path = item.get("path", "")
            if not path.lower().endswith((".tif", ".tiff")):
                continue
            blob = item["sha"]
            size = int(item.get("size", 0))
            category, include = classify_tiff(path)
            record = {
                "repo": name,
                "default_branch": branch,
                "head_commit": head,
                "path": path,
                "git_blob_sha": blob,
                "bytes": size,
                "classification": category,
                "included_in_prior_output_corpus": bool(include),
            }
            scan_paths.append(record)
            if not include:
                if "manual review required" in category:
                    excluded_review.append({"repo": name, "path": path, "git_blob_sha": blob, "bytes": size})
                continue

            prediction_repositories.add(name)
            location = _artifact_location(name, head, path, blob, size)
            entry = by_blob.get(blob)
            if entry is None:
                cache_file = args.cache / f"{blob}.tif"
                raw = cache_file.read_bytes() if cache_file.exists() else gh_raw(
                    f"repos/{OWNER}/{name}/contents/{quote(path, safe='/')}?ref={head}"
                )
                actual_blob = git_blob_sha1(raw)
                if actual_blob != blob or len(raw) != size:
                    raise RuntimeError(
                        f"raw TIFF integrity failure for {name}/{path}: size={len(raw)}/{size}, blob={actual_blob}/{blob}"
                    )
                digest = sha256(raw)
                entry = {
                    "git_blob_sha": blob,
                    "bytes": size,
                    "sha256": digest,
                    "artifacts": [],
                    "classification": category,
                }
                entries.append(entry)
                by_blob[blob] = entry
                new_unique_blobs += 1
                if not cache_file.exists():
                    cache_file.write_bytes(raw)
            locations = entry.setdefault("artifacts", [])
            if not any(
                a.get("repo") == name and a.get("ref") == head and a.get("path") == path
                for a in locations
            ):
                locations.append(location)
                new_output_locations += 1

    if truncated:
        raise RuntimeError(f"recursive GitHub trees were truncated: {truncated}")
    if tree_count != len(repos):
        raise RuntimeError(f"only {tree_count}/{len(repos)} repository trees were read")

    # Preserve stable uniqueness and make the date-scoped scan fully reviewable.
    inventory["snapshot_date_utc"] = "2026-10-06"
    inventory["repositories_scanned"] = len(repos)
    inventory["repository_names"] = [r["name"] for r in repos]
    inventory["source_repositories_with_tiff_outputs"] = sorted(prediction_repositories)
    inventory["unique_rasters"] = sorted(entries, key=lambda e: e["git_blob_sha"])
    inventory["unique_git_blobs"] = len(entries)
    inventory["unique_pinned_repo_ref_path_locations"] = sum(len(e.get("artifacts", [])) for e in entries)
    inventory["artifact_observations_before_exact_path_dedup"] = sum(len(e.get("artifacts", [])) for e in entries)
    inventory["total_unique_blob_bytes"] = sum(int(e["bytes"]) for e in entries)
    inventory["current_head_refresh"] = {
        "snapshot_date_utc": "2026-10-06",
        "owner_public_repositories_matching_gemsdoe": len(repos),
        "recursive_trees_read": tree_count,
        "recursive_trees_truncated": 0,
        "all_heads": heads,
        "current_head_tiff_paths_scanned": len(scan_paths),
        "current_head_prediction_output_locations_added_or_confirmed": new_output_locations,
        "new_unique_prediction_blobs_added": new_unique_blobs,
        "unclassified_tiff_paths_excluded_pending_manual_review": excluded_review,
        "classification_rules": "See scripts/refresh_prior_inventory.py: explicit source/template/fixture exclusions, explicit candidate/submission/model-output includes, separate proxy-circular flag.",
        "scan_paths": scan_paths,
    }
    inventory["scope"] = (
        "Union of the 2026-10-05 pinned historical prediction-output corpus and all prediction, submission, "
        "candidate, leaderboard-anchor, model-probability, calibration, archived and quarantined output TIFFs "
        "identified in all 2026-10-06 public GEMSDOE repository default-branch heads. Source rasters, labels, "
        "sample templates, proxy truth, and fixtures are classified and excluded. Paths and commits are retained; "
        "byte blobs are deduplicated by Git SHA-1. This cannot include private/unlisted repositories or future heads."
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "out": str(args.out),
        "repositories_scanned": len(repos),
        "tiff_paths_scanned": len(scan_paths),
        "current_prediction_locations_added_or_confirmed": new_output_locations,
        "new_unique_prediction_blobs_added": new_unique_blobs,
        "unique_blobs_total": len(entries),
        "included_prediction_locations_total": inventory["unique_pinned_repo_ref_path_locations"],
        "unclassified_excluded": len(excluded_review),
        "cache": str(args.cache),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
