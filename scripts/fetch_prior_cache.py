#!/usr/bin/env python3
"""Fetch the pinned public prior-raster corpus into an ignored local cache."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import urllib.request


def git_blob_sha1(data: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()


def fetch_bytes(entry: dict, owner: str) -> bytes:
    artifact = entry["artifacts"][0]
    endpoint = f"repos/{owner}/{artifact['repo']}/contents/{artifact['path']}?ref={artifact['ref']}"
    if shutil.which("gh"):
        try:
            return subprocess.check_output(["gh", "api", "-H", "Accept: application/vnd.github.raw", endpoint])
        except subprocess.CalledProcessError:
            pass
    url = artifact.get("raw_url")
    if not url:
        from urllib.parse import quote
        url = f"https://raw.githubusercontent.com/{owner}/{artifact['repo']}/{artifact['ref']}/{quote(artifact['path'], safe='/')}"
    request = urllib.request.Request(url, headers={"User-Agent": "GEMSDOE40-prior-audit/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory", type=Path, default=Path("docs/data/prior_raster_inventory.json"))
    parser.add_argument("--cache", type=Path, default=Path("/tmp/gemsdoe40-prior-cache"))
    parser.add_argument("--owner", default="buffedlizard55-lab")
    args = parser.parse_args()
    inventory = json.loads(args.inventory.read_text(encoding="utf-8"))
    args.cache.mkdir(parents=True, exist_ok=True)
    ok = 0
    for i, entry in enumerate(inventory["unique_rasters"], start=1):
        blob = entry["git_blob_sha"]
        destination = args.cache / f"{blob}.tif"
        data = destination.read_bytes() if destination.exists() else fetch_bytes(entry, args.owner)
        expected_len = int(entry["bytes"])
        actual_blob = git_blob_sha1(data)
        actual_sha256 = hashlib.sha256(data).hexdigest()
        if len(data) != expected_len or actual_blob != blob or actual_sha256 != entry["sha256"]:
            raise RuntimeError(f"integrity failure for {blob}: size={len(data)} blob={actual_blob} sha256={actual_sha256}")
        if not destination.exists():
            destination.write_bytes(data)
        ok += 1
        if i % 25 == 0 or i == len(inventory["unique_rasters"]):
            print(f"verified {i}/{len(inventory['unique_rasters'])} prior rasters")
    print(f"PASS: {ok} unique TIFF blobs cached at {args.cache}")


if __name__ == "__main__":
    main()
