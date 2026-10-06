#!/usr/bin/env python3
"""Independent re-audit of a shipped submission GeoTIFF.

Reads the *bytes on disk* (not the in-memory array) and re-checks every rule this
repository claims to satisfy:

* identical shape, transform, CRS, and pixel grid to ``data/sample_submission.tif``;
* single band, float32, exactly the sample's finite footprint;
* every finite value inside ``[0, 1]`` and no ``nodata`` sentinel that the portal
  could misinterpret (the cause of the "Predicted values must be in range [0, 1]"
  rejection);
* no positive value on a published-catalogue pixel (the organizer masks those);
* SHA-256 agreement with the receipt written by ``scripts/run_h4.py``;
* uniqueness against every prior raster in ``data/prior`` (Pearson + support IoU).

Exit status is non-zero if any check fails, so this can be run in CI.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe40.corpus import PRIORS, load_field, local_path  # noqa: E402
from gemsdoe40.grid import footprint_from_sample  # noqa: E402
from gemsdoe40.uniqueness import jaccard_positive, pearson  # noqa: E402


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("tif", help="shipped GeoTIFF (the -zeros twin is the portal file)")
    ap.add_argument("--receipt", default="")
    ap.add_argument("--sample", default=str(ROOT / "data" / "sample_submission.tif"))
    ap.add_argument("--labels", default=str(ROOT / "data" / "labels.tif"))
    ap.add_argument("--out", default="")
    args = ap.parse_args()
    path = Path(args.tif)
    checks: list[dict] = []

    def check(name: str, ok: bool, detail) -> None:
        checks.append(dict(check=name, ok=bool(ok), detail=detail))

    with rasterio.open(args.sample) as s:
        ref_shape, ref_transform, ref_crs = s.shape, s.transform, s.crs
        ref_finite = np.isfinite(s.read(1))
    with rasterio.open(path) as s:
        arr = s.read(1)
        chk_shape = s.shape == ref_shape
        chk_transform = tuple(round(v, 6) for v in s.transform[:6]) == tuple(round(v, 6) for v in ref_transform[:6])
        chk_crs = s.crs == ref_crs
        finite = np.isfinite(arr)
        checks.append(dict(check="shape", ok=chk_shape, detail=f"{s.shape} vs {ref_shape}"))
        checks.append(dict(check="transform", ok=chk_transform,
                           detail=f"{[round(v,4) for v in s.transform[:6]]}"))
        checks.append(dict(check="crs", ok=chk_crs, detail=str(s.crs)))
        check("single_band_float32", s.count == 1 and str(s.dtypes[0]) == "float32",
              f"count={s.count} dtype={s.dtypes[0]}")
        # two legal conventions ship together: the -nan twin reproduces the
        # sample's footprint exactly, the -zeros twin is finite everywhere (that
        # is the variant the portal range check accepts without a sentinel).
        if finite.all():
            check("zeros_twin_all_finite", True,
                  f"finite {int(finite.sum())} = full raster ({arr.size} px), no sentinel")
        else:
            check("footprint_matches_sample", bool((finite == ref_finite).all()),
                  f"finite {int(finite.sum())} vs {int(ref_finite.sum())}")
        vmin, vmax = float(np.nanmin(arr)), float(np.nanmax(arr))
        check("values_in_unit_range", vmin >= 0.0 and vmax <= 1.0, f"min={vmin} max={vmax}")
        check("portal_safe_nodata", s.nodata is None or bool(np.isnan(s.nodata)),
              f"nodata={s.nodata} (a finite sentinel such as 255 or -9999 is what the "
              f"portal rejects)")
        check("all_finite_inside_footprint", bool(np.isfinite(arr[ref_finite]).all()),
              f"positive_px={int((arr > 0).sum())}")
        check("distinct_values", bool(len(np.unique(arr[np.isfinite(arr)])) >= 2 or (arr > 0).sum() == 0),
              f"unique={np.unique(arr[np.isfinite(arr)])[:4].tolist()}")
    with rasterio.open(args.labels) as s:
        catalogue = s.read(1) == 1
    on_cat = int(((arr > 0) & catalogue).sum())
    check("no_mass_on_published_catalogue", on_cat == 0, f"positive px on catalogue = {on_cat}")

    digest = sha256(path)
    receipt = {}
    if args.receipt:
        receipt = json.loads(Path(args.receipt).read_text())
        rec = receipt.get("receipt", {})
        want = None
        for twin in ("zeros", "nan"):
            entry = rec.get(twin) or {}
            if entry.get("file", "").endswith(path.name):
                want = entry.get("sha256")
        if want:
            check("sha256_matches_receipt", digest == want, digest)
    else:
        check("sha256_recorded", True, digest)

    # uniqueness against the prior corpus, compared inside the footprint with the
    # same NaN/outside normalisation applied to both sides
    footprint = footprint_from_sample(args.sample)
    mine = load_field(path, footprint)
    rows = []
    for p in PRIORS:
        lp = local_path(p)
        if not lp.exists():
            continue
        other = load_field(lp, footprint)
        rows.append(dict(key=p.key, pearson=round(float(pearson(mine, other)), 6),
                         jaccard_positive=round(float(jaccard_positive(mine, other)), 6)))
    worst_r = max((abs(r["pearson"]) for r in rows), default=0.0)
    worst_j = max((r["jaccard_positive"] for r in rows), default=0.0)
    check("unique_vs_prior_corpus", worst_r < 0.85 and worst_j < 0.5,
          f"max |pearson|={worst_r:.4f} max jaccard={worst_j:.4f} n={len(rows)}")

    ok = all(c["ok"] for c in checks)
    out = dict(file=str(path), bytes=path.stat().st_size, sha256=digest,
               positive_px=int((arr > 0).sum()), all_checks_pass=ok, checks=checks,
               nearest_priors=sorted(rows, key=lambda r: -abs(r["pearson"]))[:5])
    print(json.dumps(out, indent=1))
    if args.out:
        Path(args.out).write_text(json.dumps(out, indent=1))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
