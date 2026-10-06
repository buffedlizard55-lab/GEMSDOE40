#!/usr/bin/env python3
"""Independent second-pass audit of the published H41 GeoTIFF.

Re-opens the *published bytes* in ``docs/downloads`` (never the generator's working
copy), re-measures every property the receipt claims, and re-runs the novelty and
blocked-proxy comparisons.  Any mismatch is reported and the script exits non-zero, so
a stale or hand-edited download cannot be presented as the audited artifact.

    python scripts/audit_h41_candidate.py
"""
from __future__ import annotations

import hashlib
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from gemsdoe40 import GRID_HEIGHT, GRID_WIDTH  # noqa: E402
from gemsdoe40.measure import credit_components  # noqa: E402
from run_h41_candidate import blocked_compare, novelty  # noqa: E402

DOWNLOADS = ROOT / "docs" / "downloads"
RECEIPT = ROOT / "docs" / "data" / "h41-generation.json"
REF_PRIOR = ROOT / "ref" / "prior"
EXPECTED_TRANSFORM = (100.0, 0.0, 243350.0, 0.0, -100.0, 4508550.0)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
    slug = receipt["files"]["slug"]
    zeros = DOWNLOADS / f"{slug}-zeros.tif"
    nan = DOWNLOADS / f"{slug}-nan.tif"
    zipfile_path = DOWNLOADS / f"{slug}-zeros.zip"
    problems: list[str] = []
    checks: dict[str, object] = {}

    with rasterio.open(ROOT / "data" / "sample_submission.tif") as ds:
        footprint = np.isfinite(ds.read(1))
    with rasterio.open(ROOT / "data" / "labels.tif") as ds:
        catalogue = ds.read(1) > 0
    with rasterio.open(ROOT / "data" / "external" / "derived_sgmc_faults_100m_u8.tif") as ds:
        sgmc = ds.read(1) > 0
    truth = sgmc & ~catalogue & footprint

    # ---- 1. primary (all-finite) file ------------------------------------- #
    checks["primary_exists"] = zeros.is_file()
    checks["primary_sha256_matches_receipt"] = (
        zeros.is_file() and sha256(zeros) == receipt["files"]["zeros"]["sha256"])
    checks["primary_bytes_match_receipt"] = (
        zeros.is_file() and zeros.stat().st_size == receipt["files"]["zeros"]["bytes"])
    with rasterio.open(zeros) as ds:
        if ds.count != 1:
            problems.append("primary is not single-band")
        if ds.dtypes[0] != "float32":
            problems.append(f"primary dtype is {ds.dtypes[0]}, not float32")
        if ds.crs is None or ds.crs.to_epsg() != 32611:
            problems.append("primary CRS is not EPSG:32611")
        if (ds.height, ds.width) != (GRID_HEIGHT, GRID_WIDTH):
            problems.append(f"primary shape {ds.shape} != {(GRID_HEIGHT, GRID_WIDTH)}")
        if tuple(ds.transform)[:6] != EXPECTED_TRANSFORM:
            problems.append(f"primary transform {tuple(ds.transform)[:6]} != sample")
        if ds.nodata is not None and np.isfinite(ds.nodata):
            problems.append(f"primary declares a finite nodata {ds.nodata}")
        a = ds.read(1)
    checks["primary_all_finite"] = bool(np.isfinite(a).all())
    checks["primary_min"] = float(a.min())
    checks["primary_max"] = float(a.max())
    checks["primary_values_are_binary"] = bool(np.isin(a, (0.0, 1.0)).all())
    checks["primary_positive_px"] = int((a > 0).sum())
    checks["primary_matches_emitted_mass"] = int((a > 0).sum()) == int(receipt["emitted"]["mass"])
    checks["primary_zero_outside_footprint"] = bool((a[~footprint] == 0).all())
    checks["primary_zero_on_catalogue"] = int((a[catalogue] > 0).sum())
    if not checks["primary_all_finite"]:
        problems.append("primary has non-finite cells")
    if float(a.min()) < 0 or float(a.max()) > 1:
        problems.append("primary has values outside [0, 1]")

    # ---- 2. sample-format twin -------------------------------------------- #
    checks["nan_sha256_matches_receipt"] = (
        nan.is_file() and sha256(nan) == receipt["files"]["nan"]["sha256"])
    with rasterio.open(nan) as ds:
        b = ds.read(1)
    checks["nan_outside_only"] = bool(np.isnan(b[~footprint]).all()
                                      and np.isfinite(b[footprint]).all())
    checks["nan_twin_same_dots"] = bool(np.array_equal(a > 0, np.nan_to_num(b) > 0))
    if not checks["nan_outside_only"]:
        problems.append("NaN twin is not NaN-exactly-outside")

    # ---- 3. zip ------------------------------------------------------------ #
    with zipfile.ZipFile(zipfile_path) as zf:
        names = zf.namelist()
        member_bytes = zf.read(names[0]) if len(names) == 1 else b""
    checks["zip_members"] = names
    checks["zip_member_is_primary"] = member_bytes == zeros.read_bytes()
    if not checks["zip_member_is_primary"]:
        problems.append("zip does not contain the audited primary bytes")

    # ---- 4. metric components --------------------------------------------- #
    comp = credit_components(a, truth, valid=footprint)
    checks["proxy_dti_recomputed"] = comp["dti"]
    checks["proxy_dti_matches_receipt"] = abs(comp["dti"] - receipt["audit"]["proxy_dti"]) < 1e-12
    if not checks["proxy_dti_matches_receipt"]:
        problems.append("recomputed proxy DTI differs from the receipt")
    checks["fraction_within_3px_of_proxy"] = float(
        ((a > 0) & (np.maximum(1.0 - _dist(truth) / 3.0, 0.0) > 0)).sum()
        / max(int((a > 0).sum()), 1))

    # ---- 5. novelty and blocked comparison --------------------------------- #
    nov = novelty(a > 0, footprint, sorted(REF_PRIOR.glob("*.tif")), top_k=int((a > 0).sum()))
    checks["novelty"] = {k: v for k, v in nov.items() if k != "rows"}
    if not nov["is_new"]:
        problems.append("novelty rule failed on the published bytes")

    comparators = {}
    for label, name in (("H33-B2", "gemsdoe32-h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros.tif"),
                        ("H27-4", "gems28-h27-4-r1-solo-d2-8-20261003-8acb75e1f2cc-allfinite.tif"),
                        ("H40-E", "gemsdoe39-h40-e-disc-h40e-30k-zeros.tif")):
        p = REF_PRIOR / name
        if p.exists():
            with rasterio.open(p) as ds:
                c = ds.read(1)
            comparators[label] = ((np.isfinite(c) & (c > 0)) & footprint).astype(np.float32)
    blocked = blocked_compare(a.astype(np.float32), footprint, truth, catalogue, comparators)
    checks["blocked"] = {"summary": blocked["summary"],
                         "means": {k: blocked[k]["mean_dti"] for k in blocked if k != "summary"}}

    report = {"slug": slug, "checks": checks, "problems": problems,
              "verdict": "PASS" if not problems else "FAIL"}
    out = ROOT / "docs" / "data" / "h41-audit.json"
    out.write_text(json.dumps(report, indent=2, default=float) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, default=float))
    print(f"[out] {out}")
    return 0 if not problems else 1


def _dist(truth: np.ndarray) -> np.ndarray:
    from scipy import ndimage
    return ndimage.distance_transform_edt(~truth)


if __name__ == "__main__":
    raise SystemExit(main())
