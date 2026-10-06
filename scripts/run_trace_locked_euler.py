#!/usr/bin/env python3
"""Generate the frozen H8 trace-locked Euler depth-consensus candidate.

Reads only: the hash-pinned H4 solution cloud CSV, the sample template, and exact known labels.
No proxy truth, holdout label, prior prediction or leaderboard value is opened.  Writes to the
ignored ``work/`` directory; publication happens only after the independent audit runner.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from importlib.metadata import version
import json
from pathlib import Path
import sys
import time

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from gemsdoe40.raster import validate_candidate, write_candidate, write_json  # noqa: E402
from gemsdoe40.trace_lock import CLOUD_SHA256, build_trace_locked_field, load_cloud_csv  # noqa: E402
from gemsdoe40.research_uniqueness import file_sha256  # noqa: E402
from acquire_data import PINS  # noqa: E402

PREREGISTRATION = "docs/research/h8-preregistration-20261006.md"
CLOUD_PATH = "docs/downloads/h4-euler-solutions.csv.gz"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data")
    parser.add_argument("--output", type=Path, default=ROOT / "work/h8")
    args = parser.parse_args()
    args.data = args.data.resolve()
    args.output = args.output.resolve()
    if not args.output.is_relative_to(ROOT / "work"):
        parser.error("generation output must stay under the repository's ignored work/ directory")
    started = time.monotonic()
    args.output.mkdir(parents=True, exist_ok=True)

    prereg_sha = file_sha256(ROOT / PREREGISTRATION)
    cloud_sha = file_sha256(ROOT / CLOUD_PATH)
    if cloud_sha != CLOUD_SHA256:
        raise RuntimeError(f"frozen H4 cloud hash mismatch: {cloud_sha} != {CLOUD_SHA256}")

    manifest = {
        "experiment": "H8 trace-locked Euler depth-consensus emission on the frozen H4 cloud",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "software": {"python": sys.version.split()[0],
                     **{p: version(p) for p in ("numpy", "scipy", "rasterio", "affine")}},
        "compute": "CPU; no learned weights, no GPU",
        "construction_reads_proxy_or_prior_predictions": False,
        "preregistration": {"path": PREREGISTRATION, "sha256": prereg_sha},
        "code_sha256": {str(p.relative_to(ROOT)): file_sha256(p) for p in
                        (Path(__file__), ROOT / "src/gemsdoe40/trace_lock.py",
                         ROOT / "src/gemsdoe40/contact_euler.py", ROOT / "src/gemsdoe40/raster.py")},
        "inputs": {},
    }
    for name in ("sample_submission.tif", "labels.tif"):
        path = args.data / name
        sha = file_sha256(path)
        if sha != PINS[name]:
            raise RuntimeError(f"pinned input hash mismatch: {name}")
        manifest["inputs"][name] = {"sha256": sha, "bytes": path.stat().st_size}
    manifest["inputs"]["h4-euler-solutions.csv.gz"] = {"sha256": cloud_sha,
                                                       "bytes": (ROOT / CLOUD_PATH).stat().st_size,
                                                       "origin": "frozen H4 solver output; SI=0 + offset A; "
                                                                 "docs/research/h4-preregistration-20261006.md"}

    sample = args.data / "sample_submission.tif"
    with rasterio.open(sample) as ds:
        footprint = np.isfinite(ds.read(1))
        shape, transform = ds.shape, ds.transform
    with rasterio.open(args.data / "labels.tif") as ds:
        if ds.shape != shape or ds.transform != transform:
            raise RuntimeError("labels do not match the sample grid")
        known = (ds.read(1) == 1) & footprint
    manifest["shape"] = list(shape)
    manifest["footprint_pixels"] = int(footprint.sum())
    manifest["exact_known_mask_pixels"] = int(known.sum())

    cloud = load_cloud_csv(ROOT / CLOUD_PATH, transform)
    manifest["weighted_cloud_records"] = int(len(cloud))
    field, stats = build_trace_locked_field(cloud, footprint, known)
    manifest["field_stats"] = stats

    stamp = datetime.now(timezone.utc).strftime("%Y%m%d")
    temp = args.output / "h8-candidate.tif"
    write_candidate(temp, field, sample,
                    description="GEMSDOE40 H8 trace-locked Euler depth-consensus field; research candidate")
    receipt = validate_candidate(temp, sample)
    if not receipt["valid"]:
        raise RuntimeError(f"format audit failed: {receipt['issues']}")
    filename = f"gemsdoe40-h8-tracelock-depthkde-{stamp}-{receipt['canonical_pixels_sha256'][:12]}.tif"
    destination = args.output / filename
    temp.replace(destination)
    receipt["path"] = str(destination.relative_to(ROOT))
    manifest["candidate_format"] = receipt
    manifest["positive_pixels"] = int(np.count_nonzero(field[footprint]))
    manifest["elapsed_s"] = round(time.monotonic() - started, 3)
    manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()

    np.save(args.output / "h8-field.npy", field)
    write_json(args.output / "generation.json", manifest)
    print(json.dumps({"status": "generated", "path": receipt["path"],
                      "sha256": receipt["sha256"], "positive_pixels": manifest["positive_pixels"],
                      "field_stats": stats}, indent=1))


if __name__ == "__main__":
    main()
