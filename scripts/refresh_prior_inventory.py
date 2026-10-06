#!/usr/bin/env python3
"""Conservative, pinned public-output census. Preserve history; inspect ZIPs too.

This is an owner-repository snapshot, not a claim about private/deleted artifacts.
Requires gh. All caches stay outside Git under data/prior and work/.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
from urllib.parse import quote
import zipfile

ROOT = Path(__file__).resolve().parents[1]
OWNER = "buffedlizard55-lab"


def api(endpoint: str, raw: bool = False) -> bytes | dict | list:
    command = ["gh", "api"]
    if raw:
        command += ["-H", "Accept: application/vnd.github.raw"]
    payload = subprocess.check_output(command + [endpoint], timeout=300)
    return payload if raw else json.loads(payload)


def output_reason(path: str) -> str | None:
    low = path.lower()
    if not low.endswith((".tif", ".tiff", ".zip")):
        return None
    if any(token in low for token in ("/downloads/", "downloads/", "/submissions/", "submissions/", "docs/archive/",
                                     "/scored/", "/lb_anchors/", "/legacy_candidates/", "/quarantine/",
                                     "/leaderboard_anchor/", "/suture/", "/union_po_loo/", "inputs/calibration/")):
        return "output-directory"
    base = Path(low).name
    if (low.startswith("docs/gems") or "context_detector_prob" in base or base == "prob_raw.tif"
            or (low.startswith("inputs/") and base.startswith(("gems16-", "gems19-")))):
        return "raw-model-output-or-root-docs-download"
    if any(token in base for token in ("submission", "prediction", "candidate")) and not any(
            token in base for token in ("sample_submission", "example_submission", "submission_format")):
        return "prediction-filename"
    return None


def blob_hash(data: bytes) -> str:
    return hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=ROOT / "docs/data/prior_raster_inventory.json")
    parser.add_argument("--output", type=Path, default=ROOT / "docs/data/prior-inventory-20261006.json")
    parser.add_argument("--cache", type=Path, default=ROOT / "data/prior")
    args = parser.parse_args()
    old = json.loads(args.base.read_text())
    entries = {e["git_blob_sha"]: e for e in old["unique_rasters"]}
    args.cache.mkdir(parents=True, exist_ok=True)
    repositories = []
    page = 1
    while True:
        batch = api(f"users/{OWNER}/repos?per_page=100&page={page}")
        repositories.extend(r for r in batch if re.search(r"GEMS.*DOE", r["name"], re.I))
        if len(batch) < 100:
            break
        page += 1
    observed, excluded, snapshots, errors, archives = [], [], [], [], []

    def scan(repo: dict) -> tuple:
        name = repo["name"]
        commit = api(f"repos/{OWNER}/{name}/commits/{repo['default_branch']}")["sha"]
        tree = api(f"repos/{OWNER}/{name}/git/trees/{commit}?recursive=1")
        if tree.get("truncated"):
            raise RuntimeError(f"truncated tree: {name}")
        selected, skipped = [], []
        for item in tree["tree"]:
            if item["type"] != "blob" or not item["path"].lower().endswith((".tif", ".tiff", ".zip")):
                continue
            artifact = {"repo": name, "ref": commit, "path": item["path"], "bytes": item["size"],
                        "git_blob_sha": item["sha"], "kind": output_reason(item["path"]),
                        "raw_url": f"https://raw.githubusercontent.com/{OWNER}/{name}/{commit}/{quote(item['path'], safe='/')}"}
            (selected if artifact["kind"] else skipped).append(artifact)
        return {"repo": name, "commit": commit, "tree": tree["sha"], "output_paths": len(selected)}, selected, skipped

    with ThreadPoolExecutor(max_workers=5) as pool:
        pending = {pool.submit(scan, r): r["name"] for r in repositories}
        for future in as_completed(pending):
            try:
                snap, paths, skip = future.result()
                snapshots.append(snap); observed.extend(paths); excluded.extend(skip)
                print(f"census {snap['repo']}: {len(paths)} output TIFF/ZIP paths", flush=True)
            except Exception as exc:
                errors.append({"repo": pending[future], "error": str(exc)})
    if errors:
        raise RuntimeError(f"incomplete repository census: {errors}")
    census_path = ROOT / "work/source-audit/current-census.json"
    census_path.parent.mkdir(parents=True, exist_ok=True)
    census_path.write_text(json.dumps({"observed": observed, "excluded": excluded, "snapshots": snapshots}, indent=2))

    def download(artifact: dict) -> tuple[dict, bytes]:
        blob = artifact["git_blob_sha"]
        cached = args.cache / f"{blob}.tif"
        if cached.exists() and not artifact["path"].lower().endswith(".zip"):
            data = cached.read_bytes()
        else:
            data = api(f"repos/{OWNER}/{artifact['repo']}/contents/{quote(artifact['path'], safe='/')}?ref={artifact['ref']}", raw=True)
        if len(data) != artifact["bytes"] or blob_hash(data) != blob:
            raise RuntimeError(f"size/Git blob integrity failure: {artifact}")
        return artifact, data

    with ThreadPoolExecutor(max_workers=4) as pool:
        pending = {pool.submit(download, a): a for a in observed}
        for future in as_completed(pending):
            artifact, data = future.result()
            if artifact["path"].lower().endswith(".zip"):
                with zipfile.ZipFile(io.BytesIO(data)) as z:
                    members = [p for p in z.infolist() if p.filename.lower().endswith((".tif", ".tiff"))]
                    if not members:
                        archives.append({"artifact": artifact, "sha256": hashlib.sha256(data).hexdigest(),
                                         "tiff_members": [], "status": "INVALID advertised submission archive: no TIFF",
                                         "other_members": z.namelist()})
                        print(f"IRREGULARITY: no TIFF in {artifact['repo']}/{artifact['path']}", flush=True)
                        continue
                    # Never extract paths; reject zip bombs / non-submission data bundles.
                    if any(m.file_size > 160_000_000 for m in members):
                        raise RuntimeError(f"unexpected oversized TIFF in {artifact['path']}")
                    payloads = [(z.read(m), {**artifact, "container_path": artifact["path"],
                                          "member": m.filename, "kind": "TIFF-inside-output-ZIP"}) for m in members]
                archives.append({"artifact": artifact, "sha256": hashlib.sha256(data).hexdigest(),
                                 "tiff_members": [m.filename for m in members]})
            else:
                payloads = [(data, artifact)]
            for pixels, location in payloads:
                blob = blob_hash(pixels)
                sha256 = hashlib.sha256(pixels).hexdigest()
                loc = {k: v for k, v in location.items() if k != "git_blob_sha"}
                if blob not in entries:
                    entries[blob] = {"git_blob_sha": blob, "sha256": sha256, "bytes": len(pixels), "artifacts": []}
                e = entries[blob]
                if e["sha256"] != sha256 or e["bytes"] != len(pixels):
                    raise RuntimeError(f"inconsistent prior entry: {blob}")
                if loc not in e["artifacts"]:
                    e["artifacts"].append(loc)
                dest = args.cache / f"{blob}.tif"
                if not dest.exists():
                    temp = dest.with_suffix(".partial")
                    temp.write_bytes(pixels); temp.replace(dest)
    report = {
        "snapshot_utc": datetime.now(timezone.utc).isoformat(), "owner": OWNER,
        "scope": "All current public owner GEMS*DOE repository output TIFFs and TIFFs inside output ZIPs, plus all historical pins in the 20261005 inventory. No private/deleted artifacts or later changes are claimed.",
        "parent_inventory_sha256": hashlib.sha256(args.base.read_bytes()).hexdigest(),
        "repository_snapshots": sorted(snapshots, key=lambda s: s["repo"]),
        "output_paths_observed": len(observed), "zip_archives_checked": archives,
        "excluded_source_or_unclassified_rasters": excluded,
        "unique_git_blobs": len(entries), "total_unique_blob_bytes": sum(e["bytes"] for e in entries.values()),
        "unique_rasters": sorted(entries.values(), key=lambda e: e["git_blob_sha"]),
        "errors": errors,
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS: {len(entries)} unique raw TIFF blobs; {len(archives)} ZIPs inspected; exclusions {len(excluded)} require review", flush=True)


if __name__ == "__main__":
    main()
