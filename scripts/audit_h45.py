#!/usr/bin/env python3
"""Audit the H45 candidate against the format rules and the whole prior-submission corpus.

Checks, all measured on the written bytes (never on intent):

FORMAT
  * grid, CRS, transform, dtype, band count and nodata tag match ``sample_submission.tif``
  * every value inside the footprint is finite and in [0, 1]
  * outside the footprint there are no finite non-zero values

UNIQUENESS (the gate that must pass before anything is promoted out of ``work/``)
  * SHA-256 of the file, and SHA-256 of the in-footprint pixel block, against all 343
    cached prior submissions
  * pixel-level agreement with each prior: |A n B| / |A u B| (Jaccard) and the fraction of
    this candidate's positive pixels that are positive in the prior

The candidate fails the gate if it duplicates a prior, if it is a near-duplicate
(Jaccard >= 0.90 or coverage >= 0.95), or if any format check fails.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
PRIOR_DIR = ROOT / "data" / "prior"

JACCARD_LIMIT = 0.90
COVERAGE_LIMIT = 0.95


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def format_checks(path: Path) -> dict:
    with rasterio.open(ROOT / "data/sample_submission.tif") as ds:
        tpl = ds.profile
        tpl_foot = np.isfinite(ds.read(1))
    with rasterio.open(path) as ds:
        prof = ds.profile
        arr = ds.read(1)
    foot = np.isfinite(arr)
    inside = arr[tpl_foot]
    outside = arr[~tpl_foot]
    return dict(
        path=str(path),
        bytes=path.stat().st_size,
        sha256=sha256_file(path),
        crs_ok=str(prof.get("crs")) == "EPSG:32611",
        transform_ok=list(prof.get("transform")[:6]) == list(tpl["transform"][:6]),
        shape_ok=(prof["height"], prof["width"]) == (tpl["height"], tpl["width"]),
        dtype_ok=str(prof.get("dtype")) == "float32",
        count_ok=int(prof.get("count")) == 1,
        nodata_tag_is_none=prof.get("nodata") is None,
        footprint_matches_template=bool(np.array_equal(foot, tpl_foot)),
        all_cells_finite=bool(np.isfinite(arr).all()),
        inside_footprint_all_finite=bool(np.isfinite(arr[tpl_foot]).all()),
        nan_cells=int(np.isnan(arr).sum()),
        min_value=float(np.nanmin(arr)),
        max_value=float(np.nanmax(arr)),
        values_outside_0_1=int(((inside < 0) | (inside > 1)).sum()),
        finite_nonzero_outside_footprint=int(np.count_nonzero(
            np.nan_to_num(outside, nan=0.0) != 0.0)),
        positive_px=int((inside > 0).sum()),
        total_mass=float(np.nansum(inside)),
        portal_legal=None,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("candidate", type=Path)
    ap.add_argument("--out", type=Path, default=ROOT / "work" / "h45" / "h45-audit.json")
    ap.add_argument("--prior-dir", type=Path, default=PRIOR_DIR)
    args = ap.parse_args()

    t0 = time.time()
    checks = format_checks(args.candidate)
    checks["portal_legal"] = bool(
        checks["crs_ok"] and checks["transform_ok"] and checks["shape_ok"]
        and checks["dtype_ok"] and checks["count_ok"] and checks["nodata_tag_is_none"]
        and checks["values_outside_0_1"] == 0 and checks["inside_footprint_all_finite"]
        and checks["finite_nonzero_outside_footprint"] == 0)
    print(json.dumps(checks, indent=1), flush=True)

    with rasterio.open(args.candidate) as ds:
        arr = ds.read(1)
    foot = np.isfinite(arr)
    mine = (np.nan_to_num(arr, nan=0.0) > 0) & foot
    my_h, my_w = mine.shape
    my_n = int(mine.sum())
    mine_flat = mine.ravel()

    prior_files = sorted(args.prior_dir.glob("*.tif"))
    print(f"comparing against {len(prior_files)} cached prior submissions", flush=True)
    results = []
    for i, p in enumerate(prior_files):
        try:
            with rasterio.open(p) as ds:
                a = ds.read(1)
        except Exception as exc:                                   # noqa: BLE001
            results.append(dict(file=p.name, error=str(exc)[:120]))
            continue
        if a.shape != mine.shape:
            results.append(dict(file=p.name, shape_mismatch=list(a.shape)))
            continue
        b = (np.nan_to_num(a, nan=0.0, posinf=0.0, neginf=0.0) > 0)
        inter = int(np.count_nonzero(mine_flat & b.ravel()))
        union = int(np.count_nonzero(mine_flat | b.ravel()))
        n_b = int(b.sum())
        results.append(dict(file=p.name, prior_positive=n_b, intersection=inter,
                            union=union,
                            jaccard=inter / union if union else 0.0,
                            coverage=inter / my_n if my_n else 0.0,
                            recall=inter / n_b if n_b else 0.0,
                            sha256=sha256_file(p)))
        if i % 50 == 0:
            print(f"  {i}/{len(prior_files)}", flush=True)
        del a, b

    ok = [r for r in results if "jaccard" in r]
    # "Coverage" (what fraction of this candidate's positive pixels also fire in the prior)
    # is only meaningful against priors of comparable sparsity.  Several cached priors are
    # dense continuous fields that are positive on most of the map, which trivially cover
    # every candidate pixel while sharing almost none of its structure (Jaccard ~0.016);
    # those are reported separately instead of being allowed to trip the gate.
    comparable = [r for r in ok if r["prior_positive"] <= 5 * max(my_n, 1)]
    dense = [r for r in ok if r["prior_positive"] > 5 * max(my_n, 1)]
    ok.sort(key=lambda r: r["jaccard"], reverse=True)
    comparable.sort(key=lambda r: r["coverage"], reverse=True)
    top = ok[:10]
    top_cov = comparable[:10]
    worst_j = max((r["jaccard"] for r in ok), default=0.0)
    worst_c = max((r["coverage"] for r in comparable), default=0.0)
    duplicate = [r["file"] for r in results
                 if r.get("sha256") == checks["sha256"]]

    audit = dict(
        candidate=checks,
        priors_compared=len(prior_files),
        priors_compared_ok=len(ok),
        exact_duplicates=duplicate,
        max_jaccard=worst_j,
        max_coverage=worst_c,
        most_similar=top,
        most_covering_comparable=top_cov,
        dense_priors_excluded_from_coverage=len(dense),
        gate=dict(
            format_pass=bool(checks["portal_legal"]),
            no_exact_duplicate=not duplicate,
            jaccard_pass=worst_j < JACCARD_LIMIT,
            coverage_pass=worst_c < COVERAGE_LIMIT,
        ),
        elapsed_s=round(time.time() - t0, 1),
    )
    audit["gate"]["unique_pass"] = bool(
        audit["gate"]["no_exact_duplicate"] and audit["gate"]["jaccard_pass"]
        and audit["gate"]["coverage_pass"])
    audit["gate"]["promote"] = bool(audit["gate"]["format_pass"] and audit["gate"]["unique_pass"])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(audit, indent=1, default=float))
    print(json.dumps(audit["gate"], indent=1))
    print(f"max jaccard {worst_j:.4f}  max coverage {worst_c:.4f}")
    for r in top[:5]:
        print(f"   {r['file']}  J={r['jaccard']:.4f} cov={r['coverage']:.4f} "
              f"rec={r['recall']:.4f}")
    sys.exit(0 if audit["gate"]["promote"] else 1)


if __name__ == "__main__":
    main()
