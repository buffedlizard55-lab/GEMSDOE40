#!/usr/bin/env python3
"""Uniqueness + format audit of the (withdrawn) H4-A artifact.

The full H4-A pipeline was re-run with the registered uniqueness audit enabled,
but the audit's corpus pass did not complete inside the sandbox budget, so this
script performs the same registered checks directly on the frozen artifact:

  * grid/range/nodata re-read (the same checks the submission gate applies),
  * Pearson correlation and positive-pixel Jaccard against every cached prior
    raster, using the thresholds registered in ``src/gemsdoe40/uniqueness.py``
    (near-duplicate if |Pearson| >= 0.85 or Jaccard >= 0.50).

It is an *audit*, not a promotion: the H4-A arm is closed as a negative result
(see docs/reports/h4a-negative-result-20261006.md).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe40.grid import check_submission, sha256  # noqa: E402

NEAR_DUPLICATE_CORR = 0.85
NEAR_DUPLICATE_JACCARD = 0.50


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--artifact", type=Path,
                    default=ROOT / "work" / "withdrawn" /
                    "gems40-h4a-euler-sicontact-depthcluster-20261006T012528Z-570e6300-nan.tif")
    ap.add_argument("--template", type=Path, default=ROOT / "data" / "sample_submission.tif")
    ap.add_argument("--prior-cache", type=Path, default=ROOT / "prior_cache")
    ap.add_argument("--out", type=Path, default=ROOT / "docs" / "data" / "h4a_uniqueness_audit.json")
    args = ap.parse_args()

    t0 = time.time()
    with rasterio.open(args.template) as _t:
        _foot = np.isfinite(_t.read(1))
    fmt = check_submission(args.artifact, _foot)
    with rasterio.open(args.artifact) as src:
        cand = src.read(1).astype(np.float32)
    cand_binary = cand > 0
    cand_flat = cand.ravel()
    cand_flat = cand_flat - cand_flat.mean()
    cand_norm = float(np.linalg.norm(cand_flat))
    n_cand = int(cand_binary.sum())

    files = sorted(args.prior_cache.glob("*.tif"))
    rows = []
    worst_corr, worst_jac, worst_name = 0.0, 0.0, None
    for i, path in enumerate(files, start=1):
        try:
            with rasterio.open(path) as src:
                if src.shape != cand.shape:
                    rows.append({"file": path.name, "error": f"shape {src.shape}"})
                    continue
                prior = src.read(1).astype(np.float32)
        except Exception as exc:
            rows.append({"file": path.name, "error": str(exc)})
            continue
        prior_flat = np.where(np.isfinite(prior), prior, 0.0).ravel()
        centered = prior_flat - prior_flat.mean()
        denom = cand_norm * float(np.linalg.norm(centered))
        corr = float(np.dot(cand_flat, centered) / denom) if denom else 0.0
        prior_pos = prior_flat > 0
        inter = int((cand_binary.ravel() & prior_pos).sum())
        union = int((cand_binary.ravel() | prior_pos).sum())
        jac = float(inter / union) if union else 0.0
        rows.append({"file": path.name, "pearson": round(corr, 6),
                     "jaccard_positive": round(jac, 6),
                     "prior_positive_px": int(prior_pos.sum())})
        if abs(corr) >= abs(worst_corr):
            worst_corr, worst_name = corr, path.name
        if jac >= worst_jac:
            worst_jac = jac
        if i % 50 == 0:
            print(f"{i}/{len(files)} compared ({time.time()-t0:.0f}s)", flush=True)

    is_new = (abs(worst_corr) < NEAR_DUPLICATE_CORR) and (worst_jac < NEAR_DUPLICATE_JACCARD)
    payload = {
        "arm": "H4-A",
        "status": "withdrawn — negative result on the registered comparator",
        "artifact": str(args.artifact.relative_to(ROOT)),
        "artifact_sha256": sha256(args.artifact),
        "candidate_positive_px": n_cand,
        "format_check": fmt,
        "uniqueness": {
            "n_priors": len(files),
            "n_compared": len([r for r in rows if "pearson" in r]),
            "is_new": bool(is_new),
            "worst_pearson": worst_corr,
            "worst_jaccard": worst_jac,
            "worst_file": worst_name,
            "threshold_pearson": NEAR_DUPLICATE_CORR,
            "threshold_jaccard": NEAR_DUPLICATE_JACCARD,
            "rows": rows,
        },
        "elapsed_s": round(time.time() - t0, 1),
    }
    args.out.write_text(json.dumps(payload, indent=2, default=str) + "\n")
    print(f"format ok={fmt.get('ok') if isinstance(fmt, dict) else fmt}")
    print(f"uniqueness is_new={is_new} worst |pearson|={abs(worst_corr):.4f} "
          f"(threshold {NEAR_DUPLICATE_CORR}) worst jaccard={worst_jac:.4f} "
          f"(threshold {NEAR_DUPLICATE_JACCARD}) over {payload['uniqueness']['n_compared']} priors")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
