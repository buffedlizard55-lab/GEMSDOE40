#!/usr/bin/env python3
"""Generate the frozen H13 candidate: binary emission of H8 depth-consensus crests.

Same physics as H8 stages 1-6 (frozen H4 cloud, per-family KDE, consensus field, orientation,
orientation-NMS crests, along-strike depth consensus) with exactly the H8 constants.  The only
difference is emission: 1.0 on passing crest cells, 0.0 elsewhere, NaN outside the footprint,
exact known pixels zeroed.  No smoothing, dilation, or top-k budget.
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
from gemsdoe40.trace_lock import (  # noqa: E402
    CLOUD_SHA256, consensus_field, family_kdes, load_cloud_csv,
    orientation_bins, ridge_crests, strike_consensus,
)
from gemsdoe40.research_uniqueness import file_sha256  # noqa: E402
from acquire_data import PINS  # noqa: E402

PREREGISTRATION = "docs/research/h8-preregistration-20261006.md"
CLOUD_PATH = "docs/downloads/h4-euler-solutions.csv.gz"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data")
    parser.add_argument("--output", type=Path, default=ROOT / "work/h13")
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
        "experiment": "H13 value-concentrated crest emission of the H8 depth consensus",
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
                                                       "bytes": (ROOT / CLOUD_PATH).stat().st_size}

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

    # Exact H8 stages 1-6, then the frozen H13 binary emission.
    m, g, kde_stats = family_kdes(cloud, footprint)
    consensus = consensus_field(m, g, footprint, known)
    bins, aniso = orientation_bins(consensus, footprint)
    ridge = ridge_crests(consensus, bins)
    passes, gamma, stats = strike_consensus(cloud, bins, ridge)
    if not passes.any():
        raise RuntimeError("no crest passed the along-strike depth consensus; stop, do not relax gates")
    field = np.where(passes & ~known, 1.0, 0.0).astype(np.float32)

    manifest["field_stats"] = {**kde_stats, **stats,
                               "oriented_cells": int((bins >= 0).sum()),
                               "positive_output_cells": int(np.count_nonzero(field)),
                               "known_pixels_zeroed": int(known.sum()),
                               "output_mass_sum": float(field.sum())}

    stamp = datetime.now(timezone.utc).strftime("%Y%m%d")
    temp = args.output / "h13-candidate.tif"
    write_candidate(temp, field, sample,
                    description="GEMSDOE40 H13 binary Euler depth-consensus crest emission; research candidate")
    receipt = validate_candidate(temp, sample)
    if not receipt["valid"]:
        raise RuntimeError(f"format audit failed: {receipt['issues']}")
    filename = f"gemsdoe40-h13-crest-binary-{stamp}-{receipt['canonical_pixels_sha256'][:12]}.tif"
    destination = args.output / filename
    temp.replace(destination)
    receipt["path"] = str(destination.relative_to(ROOT))
    manifest["candidate_format"] = receipt
    manifest["positive_pixels"] = int(np.count_nonzero(field))
    manifest["elapsed_s"] = round(time.monotonic() - started, 3)
    manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()

    write_json(args.output / "generation.json", manifest)
    print(json.dumps({"status": "generated", "path": receipt["path"], "sha256": receipt["sha256"],
                      "positive_pixels": manifest["positive_pixels"],
                      "field_stats": manifest["field_stats"]}, indent=1))


if __name__ == "__main__":
    main()
