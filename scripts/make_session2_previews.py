#!/usr/bin/env python3
"""Deterministic session-2 previews of the committed H13/H8 TIFFs; never AI imagery."""
from __future__ import annotations
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import rasterio
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gemsdoe40.research_uniqueness import file_sha256


def main():
    with rasterio.open(ROOT / "data/sample_submission.tif") as ds:
        foot = np.isfinite(ds.read(1))
        bounds = ds.bounds
    extent = (bounds.left / 1000, bounds.right / 1000, bounds.bottom / 1000, bounds.top / 1000)
    out = ROOT / "docs/assets"
    out.mkdir(exist_ok=True)
    cmap = LinearSegmentedColormap.from_list("gems", ["#163238", "#205956", "#409883", "#c7d998", "#f8d88e"])
    cmap.set_bad("#10282e")
    targets = [
        ("h13-raster.png", ROOT / "docs/downloads/gemsdoe40-h13-crest-binary-20261006-a5d5b80a8476.tif",
         "H13 · binary Euler depth-consensus crests"),
        ("h8-raster.png", ROOT / "docs/downloads/gemsdoe40-h8-tracelock-depthkde-20261006-373fa53b12e9.tif",
         "H8 · trace-locked Euler depth-KDE"),
    ]
    receipt = {"display": "4x4 maximum pooling for legibility; identical extent and 0-1 colors; "
                          "display-only, TIFF values unaltered", "rasters": []}
    shown_foot = ndimage.maximum_filter(foot.astype(np.uint8), size=4, mode="constant")[::4, ::4].astype(bool)
    for name, path, title in targets:
        with rasterio.open(path) as ds:
            arr = ds.read(1)
        arr = np.nan_to_num(arr, nan=0.0)
        display = ndimage.maximum_filter(arr, size=4, mode="constant")[::4, ::4]
        image = np.ma.masked_where(~shown_foot, display)
        fig, ax = plt.subplots(figsize=(6.7, 7.9), facecolor="#10282e", layout="constrained")
        ax.set_facecolor("#10282e")
        im = ax.imshow(image, origin="upper", extent=extent, cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
        ax.contour(shown_foot[::-1], [.5], extent=extent, colors=["#607d7a"], linewidths=.45)
        ax.set_title(title, color="#ebefe6", fontsize=13, loc="left", pad=14)
        ax.set_xlabel("Easting · UTM 11N (km)", color="#b3cac5", fontsize=9)
        ax.set_ylabel("Northing (km)", color="#b3cac5", fontsize=9)
        ax.tick_params(colors="#b3cac5", labelsize=8)
        for spine in ax.spines.values():
            spine.set_color("#45605f")
        ax.plot([bounds.left / 1000 + 15, bounds.left / 1000 + 65], [bounds.bottom / 1000 + 18] * 2,
                color="#d8e1d5", lw=2)
        ax.text(bounds.left / 1000 + 15, bounds.bottom / 1000 + 25, "50 km", fontsize=8, color="#d8e1d5")
        bar = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=.04, pad=.08, aspect=35)
        bar.set_label("Uncalibrated confidence · not verified fault probability", color="#b3cac5", fontsize=9)
        bar.ax.tick_params(colors="#b3cac5", labelsize=8)
        fig.savefig(out / name, dpi=150, metadata={"Description": "Measured model output, not observed faults. EPSG:32611."})
        plt.close(fig)
        receipt["rasters"].append({"file": f"docs/assets/{name}", "source_sha256": file_sha256(path),
                                   "png_sha256": file_sha256(out / name)})
    (ROOT / "docs/data/session2-previews.json").write_text(json.dumps(receipt, indent=1))
    print(json.dumps(receipt, indent=1))


if __name__ == "__main__":
    main()
