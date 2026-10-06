#!/usr/bin/env python3
"""Extend the dated prior census through reachable PUBLIC branch histories.

No branches/checkouts are created. GitHub commit/tree APIs are read-only; large
blobs remain in ignored data/prior. Public heads are pinned before traversal.
Any API/truncation/integrity failure prevents replacing the inventory.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import io
import json
import tempfile
from pathlib import Path
from urllib.parse import quote
import zipfile

from refresh_prior_inventory import ROOT, OWNER, api, output_reason, blob_hash


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh-heads", action="store_true", help="Pin then-current public branch heads for a NEW census; default resumes this dated census")
    args = parser.parse_args()
    path = ROOT / "docs/data/prior-inventory-20261006.json"
    inventory = json.loads(path.read_text())
    historical = ROOT / "work/source-audit/history"
    historical.mkdir(parents=True, exist_ok=True)
    entries = {e["git_blob_sha"]: e for e in inventory["unique_rasters"]}
    known_archives = {a["artifact"]["git_blob_sha"] for a in inventory["zip_archives_checked"]}
    pinned = inventory["repository_snapshots"]
    started = inventory.get("history_expansion", {}).get("started_utc", datetime.now(timezone.utc).isoformat())

    def cached_json(endpoint, destination):
        if destination.exists():
            return json.loads(destination.read_text())
        value = api(endpoint)
        destination.parent.mkdir(parents=True, exist_ok=True)
        # Identical Git trees can be requested concurrently from different commits.
        # Unique temporary names keep cache publication atomic without a rename race.
        with tempfile.NamedTemporaryFile(mode="w", dir=destination.parent, suffix=".partial", delete=False) as stream:
            stream.write(json.dumps(value))
            tmp = Path(stream.name)
        tmp.replace(destination)
        return value

    def history(snapshot):
        repo = snapshot["repo"]
        # Retain the dated snapshot AND every currently public branch head.
        branch_file = historical / repo / "pinned-public-heads.json"
        if branch_file.exists() and not args.refresh_heads:
            heads = json.loads(branch_file.read_text())
        else:
            branches, page = [], 1
            while True:
                batch = api(f"repos/{OWNER}/{repo}/branches?per_page=100&page={page}")
                branches.extend(batch)
                if len(batch) < 100: break
                page += 1
            heads = [{"branch": b["name"], "sha": b["commit"]["sha"]} for b in branches]
            heads.append({"branch": "dated-main-snapshot", "sha": snapshot["commit"]})
            branch_file.parent.mkdir(parents=True, exist_ok=True)
            branch_file.write_text(json.dumps(heads, indent=2))
        commits = {}
        for head in {h["sha"] for h in heads}:
            page = 1
            while True:
                batch = cached_json(f"repos/{OWNER}/{repo}/commits?sha={head}&per_page=100&page={page}",
                                    historical / repo / f"commits-{head}-{page}.json")
                for c in batch:
                    commits[c["sha"]] = c["commit"]["tree"]["sha"]
                if len(batch) < 100: break
                page += 1
        print(f"history {repo}: {len(heads)-1} public heads; {len(commits)} reachable commits", flush=True)
        return {"repo": repo, "public_heads": heads, "commits": commits}

    with ThreadPoolExecutor(max_workers=5) as pool:
        histories = list(pool.map(history, pinned))
    trees = []
    for h in histories:
        for commit, tree in h["commits"].items():
            trees.append((h["repo"], commit, tree))
    print(f"Traverse {len(trees)} commit trees, without Git checkout or branch changes", flush=True)

    def scan(item):
        repo, commit, tree_sha = item
        tree = cached_json(f"repos/{OWNER}/{repo}/git/trees/{tree_sha}?recursive=1",
                           historical / "trees" / f"{tree_sha}.json")
        if tree.get("truncated"):
            raise RuntimeError(f"truncated history tree: {repo}/{commit}")
        found = []
        for f in tree["tree"]:
            if f["type"] != "blob" or not f["path"].lower().endswith((".tif", ".tiff", ".zip")):
                continue
            found.append({"repo": repo, "ref": commit, "path": f["path"], "bytes": f["size"],
                          "git_blob_sha": f["sha"], "kind": output_reason(f["path"]),
                          "raw_url": f"https://raw.githubusercontent.com/{OWNER}/{repo}/{commit}/{quote(f['path'], safe='/')}"})
        return found

    selected, excluded = {}, {}
    with ThreadPoolExecutor(max_workers=6) as pool:
        for i, found in enumerate(pool.map(scan, trees), 1):
            for artifact in found:
                key = (artifact["repo"], artifact["path"], artifact["git_blob_sha"])
                (selected if artifact["kind"] else excluded).setdefault(key, artifact)
            if i % 100 == 0: print(f"history trees {i}/{len(trees)}", flush=True)
    discovery = {"started_utc": started, "histories": histories, "selected": list(selected.values()), "excluded": list(excluded.values())}
    (historical / "discovery.json").write_text(json.dumps(discovery, indent=2) + "\n")
    new = {}
    for a in selected.values():
        blob = a["git_blob_sha"]
        if a["path"].lower().endswith(".zip"):
            if blob not in known_archives: new.setdefault(blob, a)
        elif blob not in entries:
            new.setdefault(blob, a)
    print(f"Discovered {len(new)} uninspected historical TIFF/ZIP blobs; {len(excluded)} exclusions need review", flush=True)

    def download(artifact):
        payload = api(f"repos/{OWNER}/{artifact['repo']}/contents/{quote(artifact['path'], safe='/')}?ref={artifact['ref']}", raw=True)
        if len(payload) != artifact["bytes"] or blob_hash(payload) != artifact["git_blob_sha"]:
            raise RuntimeError(f"history download integrity failure: {artifact}")
        return artifact, payload

    baseline_path = path.with_name("prior-inventory-snapshot-20261006.json")
    baseline = json.loads(baseline_path.read_text()) if baseline_path.exists() else inventory
    old_count = len(baseline["unique_rasters"])
    old_archives = len(baseline["zip_archives_checked"])
    with ThreadPoolExecutor(max_workers=4) as pool:
        for artifact, data in pool.map(download, new.values()):
            if artifact["path"].lower().endswith(".zip"):
                with zipfile.ZipFile(io.BytesIO(data)) as z:
                    members = [m for m in z.infolist() if m.filename.lower().endswith((".tif", ".tiff"))]
                    if any(m.file_size > 160_000_000 for m in members):
                        raise RuntimeError("unexpected oversized historical output ZIP member")
                    payloads = [(z.read(m), {**artifact, "container_path": artifact["path"], "member": m.filename,
                                             "kind": "TIFF-inside-historical-output-ZIP"}) for m in members]
                    record = {"artifact": artifact, "sha256": hashlib.sha256(data).hexdigest(),
                              "tiff_members": [m.filename for m in members]}
                    if not members:
                        record.update(status="INVALID advertised output archive: no TIFF", other_members=z.namelist())
                    inventory["zip_archives_checked"].append(record)
            else:
                payloads = [(data, artifact)]
            for pixels, location in payloads:
                blob, sha = blob_hash(pixels), hashlib.sha256(pixels).hexdigest()
                loc = {k: v for k, v in location.items() if k != "git_blob_sha"}
                if blob not in entries:
                    entries[blob] = {"git_blob_sha": blob, "sha256": sha, "bytes": len(pixels), "artifacts": []}
                e = entries[blob]
                if e["sha256"] != sha or e["bytes"] != len(pixels):
                    raise RuntimeError("inconsistent historical TIFF identity")
                if loc not in e["artifacts"]: e["artifacts"].append(loc)
                dest = ROOT / "data/prior" / f"{blob}.tif"
                if dest.exists() and hashlib.sha256(dest.read_bytes()).hexdigest() != sha:
                    raise RuntimeError("corrupt existing historical cache")
                if not dest.exists():
                    tmp = dest.with_suffix(".partial"); tmp.write_bytes(pixels); tmp.replace(dest)
    # Include all discovered source references for newly recovered raw TIFFs.
    for a in selected.values():
        if a["git_blob_sha"] in entries:
            loc = {k: v for k, v in a.items() if k != "git_blob_sha"}
            if loc not in entries[a["git_blob_sha"]]["artifacts"]:
                entries[a["git_blob_sha"]]["artifacts"].append(loc)
    inventory["history_expansion"] = {
        "started_utc": started, "completed_utc": datetime.now(timezone.utc).isoformat(),
        "repositories": len(histories), "reachable_commits": len(trees),
        "public_branch_heads": sum(len(h["public_heads"])-1 for h in histories),
        "heads": [{"repo": h["repo"], "public_heads": h["public_heads"], "reachable_commits": len(h["commits"])} for h in histories],
        "new_unique_tiff_blobs": len(entries)-old_count,
        "additional_zip_blobs_inspected": len(inventory["zip_archives_checked"])-old_archives,
        "unique_output_path_blob_observations": len(selected),
        "historical_exclusions_requiring_review": list(excluded.values()),
        "completeness": "Every commit reachable from the pinned public branch heads plus the dated main snapshots was traversed. Unreachable/deleted refs, private repositories and later pushes are not knowable.",
    }
    inventory["scope"] = "Public owner GEMS*DOE output TIFFs/ZIPs across all commits reachable from the pinned public branch heads, dated main snapshots and retained older artifact pins; no claim about unreachable/private artifacts or later pushes."
    inventory["excluded_source_or_unclassified_rasters"] = [a for a in inventory["excluded_source_or_unclassified_rasters"] if output_reason(a["path"]) is None]
    inventory["unique_rasters"] = sorted(entries.values(), key=lambda e: e["git_blob_sha"])
    inventory["unique_git_blobs"] = len(entries)
    inventory["total_unique_blob_bytes"] = sum(e["bytes"] for e in entries.values())
    snapshot = path.with_name("prior-inventory-snapshot-20261006.json")
    if not snapshot.exists(): snapshot.write_bytes(path.read_bytes())
    tmp = path.with_suffix(".partial"); tmp.write_text(json.dumps(inventory, indent=2) + "\n"); tmp.replace(path)
    print(f"HISTORY COMPLETE: {len(entries)} raw TIFF blobs (+{len(entries)-old_count}), {len(trees)} commits. Review exclusions before novelty claims.", flush=True)


if __name__ == "__main__":
    main()
