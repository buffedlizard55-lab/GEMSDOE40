#!/usr/bin/env python3
"""Verify competition rasters and print a measured grid receipt."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

PINS = {
    "training_features.tif": "4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5",
    "labels.tif": "7ba308ccdc4418b31a178f4f1ef21aaa6e152e4028f2f6f64b01f7eb25ae4093",
    "sample_submission.tif": "2176d08e485aa2cd2860ce8df539db4faf4d76163b38a4dd8c30a40454d35cbc",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    data = ROOT / "data"
    report = {}
    missing = [name for name in PINS if not (data / name).is_file()]
    if missing:
        raise SystemExit(f"Missing required inputs: {missing}. Run bash scripts/download_competition_data.sh; no ready claim or pins overwritten.")
    for name, exp in PINS.items():
        p = data / name
        got = sha256(p)
        ok = got == exp
        print(f"{name}: {'OK' if ok else 'HASH MISMATCH'} {got}")
        if not ok:
            raise SystemExit(f"refuse to proceed: {name}")
        report[name] = got
    import numpy as np
    import rasterio
    with rasterio.open(data / "sample_submission.tif") as template:
        if template.count != 1 or template.crs.to_epsg() != 32611 or template.res != (100.0, 100.0):
            raise SystemExit("unexpected sample grid")
        footprint = np.isfinite(template.read(1))
        for name in ("labels.tif", "training_features.tif"):
            with rasterio.open(data / name) as ds:
                if ds.shape != template.shape or ds.transform != template.transform or ds.crs != template.crs:
                    raise SystemExit(f"unaligned input: {name}")
        print(f"aligned sample: {template.height} rows x {template.width} columns; {int(footprint.sum()):,} valid cells; {template.crs}")
    (ROOT / "evidence").mkdir(exist_ok=True)
    (ROOT / "evidence" / "data_pins.json").write_text(json.dumps(report, indent=2) + "\n")
    print("READY for CPU Euler generation; no neural-network training or GPU requirement implied")


if __name__ == "__main__":
    main()
