#!/usr/bin/env python3
"""Atomically acquire hash-pinned public mirrors; never claim official authentication."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
OWNER = "buffedlizard55-lab"
FEATURE_REF = "dcbbb192e56b2b32c0a131eba791dc363305d4a3"
SMALL_REF = "07345ea0604953d7efb858d9cfbc21e20c7aca0b"
PINS = {
    "training_features.tif": "4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5",
    "labels.tif": "7ba308ccdc4418b31a178f4f1ef21aaa6e152e4028f2f6f64b01f7eb25ae4093",
    "sample_submission.tif": "2176d08e485aa2cd2860ce8df539db4faf4d76163b38a4dd8c30a40454d35cbc",
    "external/derived_sgmc_faults_100m_u8.tif": "26d142c4c93282cd94f6950ab96f22aeff59fbbea523d43d662e76fa1b161b5c",
}
PART_BLOBS = [
    (94371840, "f34e5143d5e7c10014b5cfe01c85d01bcede4c0a"),
    (94371840, "4bd98099ae505f63cc6c3733954f259974efd290"),
    (94371840, "d6f883ba8822207ac2765445737ecd2ecdfd2914"),
    (94371840, "976977b4ac84c4c9b6fec72495c4e512341de872"),
    (41425484, "41a7bb7b7603c1e7669d8cc8405cfeb5efb5cf0c"),
]


def digest(path: Path, git_blob: bool = False) -> str:
    h = hashlib.sha1() if git_blob else hashlib.sha256()
    if git_blob:
        h.update(f"blob {path.stat().st_size}\0".encode())
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(endpoint: str, destination: Path, *, sha256: str | None = None,
          blob: str | None = None, size: int | None = None) -> None:
    def valid(path: Path) -> bool:
        return (path.is_file() and (size is None or path.stat().st_size == size)
                and (sha256 is None or digest(path) == sha256)
                and (blob is None or digest(path, True) == blob))
    if valid(destination):
        print(f"verified existing {destination.name}", flush=True)
        return
    if destination.exists():
        raise RuntimeError(f"existing input is corrupt; preserving rather than overwriting: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=destination.parent, suffix=".partial", delete=False) as out:
        temporary = Path(out.name)
        try:
            subprocess.run(["gh", "api", "-H", "Accept: application/vnd.github.raw", endpoint],
                           stdout=out, check=True, timeout=300)
            out.flush()
            if not valid(temporary):
                raise RuntimeError(f"download integrity check failed: {endpoint}")
            temporary.replace(destination)
        finally:
            temporary.unlink(missing_ok=True)
    print(f"acquired + verified {destination.name}", flush=True)


def acquire(data: Path, proxy: bool = True) -> dict:
    if not shutil.which("gh"):
        raise RuntimeError("gh is required for the public pinned mirror in this environment")
    data.mkdir(parents=True, exist_ok=True)
    target = data / "training_features.tif"
    if target.exists() and digest(target) != PINS[target.name]:
        raise RuntimeError("existing training_features.tif hash mismatch; no append or silent replacement")
    if not target.exists():
        partdir = data / "cache" / "bridge"
        def part(i: int) -> Path:
            size, blob = PART_BLOBS[i]
            name = f"gems-geodawn-numerical-features.tif.part-{i:03d}"
            out = partdir / name
            fetch(f"repos/{OWNER}/GEMSDOE/contents/data/bridge/{name}?ref={FEATURE_REF}",
                  out, blob=blob, size=size)
            return out
        with ThreadPoolExecutor(max_workers=3) as pool:
            parts = list(pool.map(part, range(len(PART_BLOBS))))
        temporary = target.with_suffix(".assembling")
        try:
            with temporary.open("wb") as out:
                for path in parts:
                    with path.open("rb") as inp:
                        shutil.copyfileobj(inp, out)
            if digest(temporary) != PINS[target.name]:
                raise RuntimeError("assembled feature hash mismatch")
            temporary.replace(target)
        finally:
            temporary.unlink(missing_ok=True)
    sources = {
        "labels.tif": "data/bridge/existing_faults.tif",
        "sample_submission.tif": "data/bridge/example_submission.tif",
    }
    for name, path in sources.items():
        fetch(f"repos/{OWNER}/GEMSDOE24/contents/{path}?ref={SMALL_REF}", data / name, sha256=PINS[name])
    if proxy:
        # Same basename in GEMSDOE24 is a DIFFERENT raster. Never substitute it.
        name = "external/derived_sgmc_faults_100m_u8.tif"
        sources[name] = "data/" + name
        fetch(f"repos/{OWNER}/GEMSDOE30/contents/data/{name}?ref=1f9ac110d5f1a4fac7957fe814267a0be964dd20",
              data / name, sha256=PINS[name])
    names = ["training_features.tif", *sources]
    receipt = {
        "provenance": "Public sibling-repository mirrors; pins establish consistency, NOT organizer-authenticated provenance.",
        "official_data_url": "https://www.drivendata.org/competitions/306/competition-doe-gems/data/",
        "feature_ref": FEATURE_REF, "template_and_labels_ref": SMALL_REF,
        "proxy_repo": "GEMSDOE30", "proxy_ref": "1f9ac110d5f1a4fac7957fe814267a0be964dd20",
        "proxy_path": "data/external/derived_sgmc_faults_100m_u8.tif",
        "proxy_irregularity": "GEMSDOE24 and 13GEMSDOE contain different rasters under this same filename; only GEMSDOE30 matches the frozen research proxy SHA-256.",
        "files": {name: {"sha256": digest(data / name), "bytes": (data / name).stat().st_size} for name in names},
        "proxy_caveat": "Owner-derived SGMC raster, not independently rebuilt USGS data or organizer truth.",
    }
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data")
    parser.add_argument("--no-proxy", action="store_true")
    args = parser.parse_args()
    receipt = acquire(args.data, proxy=not args.no_proxy)
    destination = ROOT / "docs/data/acquisition-20261006.json"
    destination.write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"PASS: {len(receipt['files'])} hash-verified inputs; {destination}")


if __name__ == "__main__":
    main()
