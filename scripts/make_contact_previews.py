#!/usr/bin/env python3
"""Deterministic visualizations of real raster/point values, never AI imagery."""
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
    release = json.loads((ROOT / "docs/data/current-candidate.json").read_text())
    with rasterio.open(ROOT / "data/sample_submission.tif") as ds:
        foot = np.isfinite(ds.read(1)); bounds = ds.bounds; transform = ds.transform
    extent = (bounds.left/1000, bounds.right/1000, bounds.bottom/1000, bounds.top/1000)
    out = ROOT / "docs/assets"; out.mkdir(exist_ok=True)
    cmap = LinearSegmentedColormap.from_list("gems", ["#163238", "#205956", "#409883", "#c7d998", "#f8d88e"])
    cmap.set_bad("#10282e")
    targets = [("h4-raster.png", ROOT / release["path"], "H4 · contact Euler depth-KDE"),
               ("h33-raster.png", ROOT / "data/prior/17a76895f68174cc93f3cb1597d25686f6d6bc67.tif", "H33-B2 · historical dotted field"),
               ("h2b-raster.png", ROOT / "docs/downloads/gemsdoe40-h2b-tmi-euler-natural-support-20261005.tif", "H2-B · historical support raster")]
    receipt = {"display": "4x4 maximum pooling for legibility; identical extent and 0–1 colors; display-only, not altered TIFF values", "rasters": []}
    for name, path, title in targets:
        with rasterio.open(path) as ds:
            arr = ds.read(1)
        arr = np.nan_to_num(arr, nan=0)
        display = ndimage.maximum_filter(arr, size=4, mode="constant")[::4, ::4]
        shown_foot = ndimage.maximum_filter(foot.astype(np.uint8), size=4, mode="constant")[::4, ::4].astype(bool)
        image = np.ma.masked_where(~shown_foot, display)
        fig, ax = plt.subplots(figsize=(6.7, 7.9), facecolor="#10282e", layout="constrained")
        ax.set_facecolor("#10282e")
        im = ax.imshow(image, origin="upper", extent=extent, cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
        ax.contour(shown_foot[::-1], [.5], extent=extent, colors=["#607d7a"], linewidths=.45)
        ax.set_title(title, color="#ebefe6", fontsize=13, loc="left", pad=14)
        ax.set_xlabel("Easting · UTM 11N (km)", color="#b3cac5", fontsize=9)
        ax.set_ylabel("Northing (km)", color="#b3cac5", fontsize=9)
        ax.tick_params(colors="#b3cac5", labelsize=8)
        for spine in ax.spines.values(): spine.set_color("#45605f")
        ax.plot([bounds.left/1000+15, bounds.left/1000+65], [bounds.bottom/1000+18]*2, color="#d8e1d5", lw=2)
        ax.text(bounds.left/1000+15, bounds.bottom/1000+25, "50 km", fontsize=8, color="#d8e1d5")
        bar = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=.04, pad=.08, aspect=35)
        bar.set_label("Uncalibrated confidence · not verified fault probability", color="#b3cac5", fontsize=9)
        bar.ax.tick_params(colors="#b3cac5", labelsize=8)
        fig.savefig(out / name, dpi=150, metadata={"Description": "Measured model output, not observed faults. EPSG:32611."})
        plt.close(fig)
        receipt["rasters"].append({"file": f"docs/assets/{name}", "source_sha256": file_sha256(path), "png_sha256": file_sha256(out/name)})
    # The actual depth-labelled solutions; height-corrected effective depths.
    fig, ax = plt.subplots(figsize=(6.7, 7.9), facecolor="#10282e", layout="constrained")
    ax.set_facecolor("#10282e")
    for family in ("tmi", "iso_grav_anom"):
        c = np.load(ROOT / f"work/h4/{family}-solutions.npz")["cloud"]
        c = c[c["weight"] > 0]
        east = (transform.c + (c["col"] + .5)*transform.a)/1000
        north = (transform.f + (c["row"] + .5)*transform.e)/1000
        im = ax.scatter(east, north, c=c["depth_m"]/1000, s=.9, cmap="viridis_r", vmin=0, vmax=5, alpha=.85, rasterized=True)
    ax.set_xlim(extent[:2]); ax.set_ylim(extent[2:]); ax.set_aspect("equal")
    ax.set_title("H4 · clustered solution depths", loc="left", color="#ebefe6", pad=14)
    ax.set_xlabel("Easting · UTM 11N (km)", color="#b3cac5", fontsize=9)
    ax.set_ylabel("Northing (km)", color="#b3cac5", fontsize=9)
    ax.tick_params(colors="#b3cac5", labelsize=8)
    for sp in ax.spines.values(): sp.set_color("#45605f")
    bar = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=.04, pad=.08, aspect=35)
    bar.set_label("Effective depth below reference plane (km) · color clipped at 5 km", color="#b3cac5", fontsize=8)
    bar.ax.tick_params(colors="#b3cac5", labelsize=8)
    fig.savefig(out / "h4-depth-cloud.png", dpi=150)
    plt.close(fig)
    (ROOT / "docs/data/h4-preview-receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    main()
