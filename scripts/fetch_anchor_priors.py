#!/usr/bin/env python3
"""Fetch the score-anchored prior rasters (and named comparators) from the pinned public corpus.

Every file is listed in ``docs/data/prior_raster_inventory.json`` with a git blob sha1,
a sha256 and a byte count; this script refuses to write a file whose bytes do not match
the inventory.  Output goes to the ignored cache ``ref/prior/``.

Usage:
    python scripts/fetch_anchor_priors.py            # anchors from ref/meta/prior_index.json
    python scripts/fetch_anchor_priors.py --name h33-h33-2-b2 ...
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
OWNER = "buffedlizard55-lab"
INVENTORY = ROOT / "docs" / "data" / "prior_raster_inventory.json"
INDEX = ROOT / "ref" / "meta" / "prior_index.json"


def git_blob_sha1(data: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()


def load_inventory() -> list[dict]:
    inv = json.loads(INVENTORY.read_text(encoding="utf-8"))
    return list(inv["unique_rasters"])


def basename(entry: dict) -> str:
    return Path(entry["artifacts"][0]["path"]).name


def matches(entry: dict, needles: list[str]) -> bool:
    name = basename(entry).lower()
    return any(n.lower() in name for n in needles)


def fetch(entry: dict, destination: Path) -> dict:
    blob = entry["git_blob_sha"]
    sha = entry["sha256"]
    size = int(entry["bytes"])
    if destination.exists():
        data = destination.read_bytes()
        if len(data) == size and git_blob_sha1(data) == blob and hashlib.sha256(data).hexdigest() == sha:
            return {"name": destination.name, "status": "cached"}
        raise RuntimeError(f"cached file fails inventory pins: {destination}")
    # A blob can be reachable from several pinned paths; the inventory's byte pins
    # describe the blob, so any of its listed locations reproduces the same bytes.
    art = entry["artifacts"][0]
    endpoint = f"repos/{OWNER}/{art['repo']}/contents/{art['path']}?ref={art['ref']}"
    payload = subprocess.check_output(
        ["gh", "api", "-H", "Accept: application/vnd.github.raw", endpoint], timeout=600)
    if len(payload) != size or git_blob_sha1(payload) != blob or hashlib.sha256(payload).hexdigest() != sha:
        raise RuntimeError(f"download fails inventory pins: {art['repo']}/{art['path']}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = destination.with_suffix(".partial")
    tmp.write_bytes(payload)
    tmp.replace(destination)
    return {"name": destination.name, "status": "fetched", "repo": art["repo"],
            "path": art["path"], "ref": art["ref"], "sha256": sha, "bytes": size}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, default=ROOT / "ref" / "prior")
    parser.add_argument("--name", action="append", default=[],
                        help="extra substring match on the prior filename")
    parser.add_argument("--receipt", type=Path,
                        default=ROOT / "ref" / "meta" / "anchor_fetch_receipt.json")
    args = parser.parse_args()

    index = json.loads(INDEX.read_text(encoding="utf-8"))
    wanted = [Path(meta["path"]).name for meta in index.values()]
    wanted += list(args.name)
    inv = load_inventory()
    done, missing = [], []
    for needle in wanted:
        hits = [e for e in inv if basename(e) == needle]
        if not hits:
            hits = [e for e in inv if matches(e, [needle])]
        if not hits:
            missing.append(needle)
            continue
        for entry in hits:
            done.append(fetch(entry, args.cache / basename(entry)))

    # Write the live-score association next to the bytes so the inversion can be audited.
    scores = {}
    for key, meta in index.items():
        name = Path(meta["path"]).name
        if (args.cache / name).exists():
            scores[name] = meta["live_score"]
    receipt = {
        "owner": OWNER,
        "inventory": str(INVENTORY.relative_to(ROOT)),
        "inventory_sha256": hashlib.sha256(INVENTORY.read_bytes()).hexdigest(),
        "files": done,
        "live_scores": scores,
        "missing_from_inventory": missing,
        "scope": ("public pinned mirrors, byte-verified against the inventory; this is a "
                  "reproducible public corpus, not an organizer-authenticated delivery"),
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(f"PASS: {len(done)} files cached at {args.cache}; receipt {args.receipt}")
    if missing:
        print(f"WARNING: not present in the inventory: {missing}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
