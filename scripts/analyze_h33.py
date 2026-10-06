#!/usr/bin/env python3
"""Measure what changed in the owner-reported best raster; no hidden-label inference."""
from __future__ import annotations
import json
from pathlib import Path
import sys

import numpy as np
import rasterio
from scipy.ndimage import distance_transform_edt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gemsdoe40.research_uniqueness import file_sha256
from gemsdoe40.raster import write_json


def main():
    inv = {e["git_blob_sha"]: e for e in json.loads((ROOT / "docs/data/prior-inventory-20261006.json").read_text())["unique_rasters"]}
    def read(blob):
        p = ROOT / "data/prior" / (blob + ".tif")
        if file_sha256(p) != inv[blob]["sha256"]:
            raise ValueError("baseline hash mismatch")
        with rasterio.open(p) as ds:
            a = ds.read(1)
        return a, {"git_blob_sha": blob, "sha256": file_sha256(p), "artifact": inv[blob]["artifacts"][0]}
    h33, h33id = read("17a76895f68174cc93f3cb1597d25686f6d6bc67")
    base, baseid = read("12b0a4ddf1cf76f2552621b8a8b1a1b88d783786")
    with rasterio.open(ROOT / "data/labels.tif") as ds:
        labels = ds.read(1) == 1
    with rasterio.open(ROOT / "data/sample_submission.tif") as ds:
        foot = np.isfinite(ds.read(1))
    distance = distance_transform_edt(~labels)
    b, h = (base > 0) & foot, (h33 > 0) & foot
    removed, added = b & ~h, h & ~b
    predicted_from_rule = np.where((base > 0) & (distance > 2) & foot, 1.0, 0.0).astype(np.float32)
    report = {
        "h33": h33id, "base": baseid, "date_utc": "2026-10-06",
        "score_attribution": "0.2778 for this file and 0.2708 for the base are owner-reported, not a TIFF-level organizer receipt.",
        "h33_positive_pixels": int(h.sum()), "base_positive_pixels": int(b.sum()),
        "removed_pixels": int(removed.sum()), "added_pixels": int(added.sum()),
        "mass_removed_fraction": float(removed.sum() / b.sum()),
        "removed_distance_to_known_m_min_max": [float(distance[removed].min()*100), float(distance[removed].max()*100)],
        "h33_min_distance_to_known_m": float(distance[h].min()*100),
        "exactly_base_pruned_at_distance_le_200m": bool(np.array_equal(predicted_from_rule[foot], h33[foot])),
        "conclusion": "Measured geometry is pruning, not discovery of additional fault positions. The score mechanism cannot be established causally without the test labels.",
        "metric_identity": "DTI=T/(0.2*T+0.2*FP+0.8*G); after removing predictions with TP loss L and FP saving S, score improves iff (1-0.2*DTI)*L < 0.2*DTI*S.",
        "tp_loss_per_fp_saved_limit_at_owner_best": .2*.2778/(1-.2*.2778),
        "hypothetical_tp_increase_to_03345_holding_fp_and_g_fixed": (.3345/(1-.2*.3345))/(.2778/(1-.2*.2778))-1,
        "not_a_score_forecast": True,
    }
    write_json(ROOT / "docs/data/h33-measured-analysis.json", report)
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()
