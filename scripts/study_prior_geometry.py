#!/usr/bin/env python3
"""Geometry and surrogate performance of every prior submission in the corpus.

Answers, with measurements rather than narrative:

* what the best-scoring prior (h33-2-b2, owner-reported 0.2778) actually contains -
  dot geometry, distance to the published catalogue, and how it performs on the
  surrogate instruments that cannot be gamed by re-drawing the catalogue;
* whether the owner-reported score ladder is explained by emitted mass, by
  catalogue overlap, or by the surrogate credits;
* whether any prior is a near-duplicate of another (support IoU, correlation).

Writes ``docs/data/prior_geometry.json``.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import rasterio
from scipy import ndimage
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe40.corpus import PRIORS, load_field, local_path  # noqa: E402
from gemsdoe40.grid import footprint_from_sample, read_labels  # noqa: E402
from gemsdoe40.instrument import (  # noqa: E402
    binary_dti_kdtree, isolated_catalogue_components, sgmc_off_catalogue,
)
from gemsdoe40.uniqueness import jaccard_positive, pearson  # noqa: E402

BANDS = (1, 2, 3, 5, 10, 20, 50, 100)


def geometry(pos: np.ndarray) -> dict:
    yy, xx = np.nonzero(pos)
    pts = np.column_stack([yy, xx]).astype(np.float64)
    if pts.shape[0] < 2:
        return dict(n_px=int(pts.shape[0]))
    tree = cKDTree(pts)
    d, _ = tree.query(pts, k=2)
    nn = d[:, 1]
    lab, ncomp = ndimage.label(pos, structure=np.ones((3, 3), bool))
    sizes = np.bincount(lab.ravel())[1:]
    return dict(
        n_px=int(pts.shape[0]),
        nn_spacing_px=dict(min=float(nn.min()), p10=float(np.percentile(nn, 10)),
                           median=float(np.median(nn)), p90=float(np.percentile(nn, 90)),
                           max=float(nn.max())),
        components=dict(n=int(ncomp), median_size=float(np.median(sizes)),
                        max_size=int(sizes.max()),
                        frac_in_components_gt_50px=float(sizes[sizes > 50].sum() / sizes.sum())),
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "docs" / "data" / "prior_geometry.json"))
    args = ap.parse_args()

    footprint = footprint_from_sample(ROOT / "data" / "sample_submission.tif")
    catalogue = read_labels(ROOT / "data" / "labels.tif")
    with rasterio.open(ROOT / "data" / "external" / "sgmc" / "derived_sgmc_faults_100m_u8.tif") as s:
        sgmc = s.read(1)
    off = sgmc_off_catalogue(sgmc, catalogue, footprint)
    iso = isolated_catalogue_components(catalogue, footprint, min_separation_px=6, min_size_px=12)
    instruments = {"sgmc_off_catalogue": dict(truth=off, mask=catalogue),
                   "isolated_components": dict(truth=iso["truth"], mask=iso["mask"])}
    dist_cat = ndimage.distance_transform_edt(~catalogue)

    report = {"generated_utc": None, "instruments": {
        k: dict(truth_px=int(v["truth"].sum())) for k, v in instruments.items()},
        "priors": []}
    fields: dict[str, np.ndarray] = {}
    for p in PRIORS:
        path = local_path(p)
        if not path.exists():
            continue
        field = load_field(path, footprint)
        fields[p.key] = field
        pos = field > 0
        with rasterio.open(path) as src:
            arr = src.read(1)
            meta = dict(dtype=str(src.dtypes[0]), nodata=src.nodata, crs=str(src.crs),
                        shape=list(src.shape), transform=list(src.transform)[:6])
        finite = np.isfinite(arr)
        row = dict(
            key=p.key, file=path.name, official=p.score,
            official_source="user_prompt_unverified", note=p.notes, emitted_px=int(pos.sum()),
            value_min=float(np.nanmin(arr)), value_max=float(np.nanmax(arr)),
            finite_frac_of_footprint=float(finite[footprint].mean()),
            nan_outside=bool((~finite & ~footprint).all()),
            values_binary=bool(np.all((arr[pos] == 1.0))),
            meta=meta, geometry=geometry(pos),
            catalogue_overlap=dict(
                px_on_catalogue=int((pos & catalogue).sum()),
                frac_within_1px=float((dist_cat[pos] <= 1).mean()) if pos.any() else None,
                frac_within_2px=float((dist_cat[pos] <= 2).mean()) if pos.any() else None,
                frac_within_3px=float((dist_cat[pos] <= 3).mean()) if pos.any() else None,
                frac_within_5px=float((dist_cat[pos] <= 5).mean()) if pos.any() else None,
                frac_within_10px=float((dist_cat[pos] <= 10).mean()) if pos.any() else None,
                frac_far_gt20px=float((dist_cat[pos] > 20).mean()) if pos.any() else None),
        )
        for iname, cfg in instruments.items():
            r = binary_dti_kdtree(pos, cfg["truth"], cfg["mask"], footprint)
            row[iname] = dict(dti=r["dti"], mean_credit=r["mean_credit"], kappa=r["kappa"],
                              tp=r["tp"], fp=r["fp"], fn=r["fn"])
        report["priors"].append(row)
        print(f"{p.key:34s} n={row['emitted_px']:7d} score={p.score:<7} "
              f"<=2px={row['catalogue_overlap']['frac_within_2px'] if pos.any() else float('nan'):.3f} "
              f"sgmc_h={row['sgmc_off_catalogue']['mean_credit']:.4f}", flush=True)

    # pairwise near-duplicate check
    keys = list(fields)
    pairs = []
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            fa, fb = fields[a], fields[b]
            if fa.shape != fb.shape:
                continue
            c = pearson(fa, fb)
            pairs.append(dict(a=a, b=b, pearson=round(float(c), 6),
                              jaccard=round(float(jaccard_positive(fa, fb)), 6)))
    pairs.sort(key=lambda r: -abs(r["pearson"]))
    report["near_duplicates"] = pairs[:15]

    # is the score ladder explained by mass or by surrogate credit?
    import math
    def spearman(x, y):
        rx, ry = np.argsort(np.argsort(x)), np.argsort(np.argsort(y))
        rx, ry = rx - rx.mean(), ry - ry.mean()
        return float((rx @ ry) / math.sqrt((rx @ rx) * (ry @ ry)))

    got = [r for r in report["priors"] if r["official"] is not None and r["official"] > 0.01
           and r["official_source"] == "user_prompt_unverified"]
    report["analysis"] = {}
    if len(got) >= 5:
        sc = np.array([r["official"] for r in got])
        report["analysis"] = dict(
            n=len(got),
            spearman_vs_emitted_px=spearman(sc, np.array([r["emitted_px"] for r in got])),
            spearman_vs_log_emitted_px=spearman(sc, np.log([r["emitted_px"] for r in got])),
            spearman_vs_sgmc_credit=spearman(sc, np.array([r["sgmc_off_catalogue"]["mean_credit"]
                                                           for r in got])),
            spearman_vs_frac_within_2px=spearman(sc, np.array([
                r["catalogue_overlap"]["frac_within_2px"] for r in got])),
            spearman_vs_iso_credit=spearman(sc, np.array([r["isolated_components"]["mean_credit"]
                                                          for r in got])),
        )
    Path(args.out).write_text(json.dumps(report, indent=1))
    print("analysis:", json.dumps(report["analysis"], indent=1))
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
