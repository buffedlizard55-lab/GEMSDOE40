#!/usr/bin/env python3
"""Verify an independent H4 rerun against the byte-identified released result."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(work: Path, destination: Path) -> dict:
    expected = json.loads((ROOT / "docs/data/current-candidate.json").read_text())
    generation = json.loads((work / "generation.json").read_text())
    candidate = ROOT / generation["candidate_format"]["path"]
    cloud = work / "h4-euler-solutions.csv.gz"
    checks = {
        "tiff_byte_identical": sha(candidate) == expected["sha256"],
        "canonical_pixels_identical": generation["candidate_format"]["canonical_pixels_sha256"] == expected["canonical_pixels_sha256"],
        "cloud_byte_identical": sha(cloud) == expected["cloud"]["sha256"],
        "format_pass": generation["candidate_format"]["valid"] is True,
        "no_prior_predictions_or_proxy_in_generation": generation["construction_uses_proxy_or_prior_predictions"] is False,
    }
    result = {"verified_utc": datetime.now(timezone.utc).isoformat(), "scope": "Independent computation, not a copy of the released TIFF/cloud; same CPU software and frozen parameters",
              "working_directory": str(work.relative_to(ROOT)), "expected_tiff_sha256": expected["sha256"],
              "actual_tiff_sha256": sha(candidate), "expected_cloud_sha256": expected["cloud"]["sha256"],
              "actual_cloud_sha256": sha(cloud), "software": generation.get("software"),
              "elapsed_s": generation["elapsed_s"], "code_sha256": generation["code_sha256"],
              "checks": checks, "pass": all(checks.values()),
              "caveat": "Bitwise reproduction is measured for this environment, not promised for every BLAS/FFT/GDAL platform. A future mismatch must be investigated, not hidden by renaming."}
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2) + "\n")
    if not result["pass"]:
        raise RuntimeError("reproduction differs from the published bytes; inspect receipt, do not claim identical")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, default=ROOT / "work/h4-reproduction")
    parser.add_argument("--output", type=Path, default=ROOT / "docs/data/h4-reproduction.json")
    args = parser.parse_args()
    result = verify(args.work.resolve(), args.output.resolve())
    print(json.dumps({"pass": result["pass"], "checks": result["checks"]}, indent=2))


if __name__ == "__main__":
    main()
