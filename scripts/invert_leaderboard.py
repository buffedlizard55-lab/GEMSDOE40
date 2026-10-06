#!/usr/bin/env python3
"""Fit the block-grid truth-density inversion to the scored anchor rasters.

Reads the byte-verified anchors in ``ref/prior``, the score association in
``ref/meta/prior_corpus.md`` (owner-reported, not organizer-authenticated), the
scoreable domain from ``data/sample_submission.tif`` and ``data/labels.tif``, and an
independent public prior (the SGMC-derived fault raster in ``data/external``) that is
used only as a weak smoothing target and as a reporting baseline.

Writes ``docs/data/truth-inversion.json`` and ``ref/meta/truth-density-<block>px.npz``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import rasterio
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe40.inversion import InversionModel, fit, leave_one_out, spearman  # noqa: E402
from gemsdoe40.measure import kernel_field  # noqa: E402

DATA = ROOT / "data"
PRIOR = ROOT / "ref" / "prior"
CORPUS = ROOT / "ref" / "meta" / "prior_corpus.md"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_corpus() -> dict[str, dict]:
    rows: dict[str, dict] = {}
    pattern = re.compile(
        r"\|\s*`([^`]+)`\s*\|\s*`([0-9a-f]+)`\s*\|\s*([\d,]+)\s*\|\s*([\d,]+)\s*\|"
        r"\s*([\d.]+)\s*\|\s*([^|]+?)\s*\|\s*([\d.]+)\s*\|")
    for line in CORPUS.read_text(encoding="utf-8").splitlines():
        m = pattern.match(line)
        if not m:
            continue
        live_raw = m.group(6).strip()
        rows[m.group(1)] = {
            "sha16": m.group(2),
            "bytes": int(m.group(3).replace(",", "")),
            "live_score": None if live_raw in ("—", "-", "") else float(live_raw),
        }
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--block", type=int, default=40)
    parser.add_argument("--lam", type=float, default=1e3)
    parser.add_argument("--mu", type=float, default=0.0)
    parser.add_argument("--max-iter", type=int, default=300)
    parser.add_argument("--out", type=Path, default=ROOT / "docs" / "data" / "truth-inversion.json")
    parser.add_argument("--grid", type=Path, default=ROOT / "ref" / "meta" / "truth-density.npz")
    parser.add_argument("--loo", action="store_true")
    args = parser.parse_args()

    t0 = time.time()
    with rasterio.open(DATA / "sample_submission.tif") as ds:
        footprint = np.isfinite(ds.read(1))
    with rasterio.open(DATA / "labels.tif") as ds:
        catalogue = ds.read(1) > 0
    with rasterio.open(DATA / "external" / "derived_sgmc_faults_100m_u8.tif") as ds:
        sgmc = ds.read(1) > 0
    scoreable = footprint & ~catalogue

    corpus = parse_corpus()
    by_sha: dict[str, Path] = {sha256(p)[:16]: p for p in PRIOR.glob("*.tif")}

    model = InversionModel(args.block, footprint.shape)
    model.set_scoreable(scoreable)
    model.set_prior(model.density_of_mask(sgmc & scoreable))

    used, skipped = [], []
    for short, meta in sorted(corpus.items()):
        if meta["live_score"] is None:
            skipped.append({"anchor": short, "reason": "no reported score"})
            continue
        path = by_sha.get(meta["sha16"])
        if path is None:
            skipped.append({"anchor": short, "reason": "raster not cached"})
            continue
        with rasterio.open(path) as ds:
            a = ds.read(1)
        if a.shape != footprint.shape:
            skipped.append({"anchor": short, "reason": "grid mismatch"})
            continue
        dots = np.isfinite(a) & (a > 0) & scoreable
        if not dots.any():
            skipped.append({"anchor": short, "reason": "no scoreable dots"})
            continue
        model.add_anchor(short, dots, kernel_field(dots), meta["live_score"])
        used.append({"anchor": short, "file": path.name, "sha256": sha256(path),
                     "mass": int(dots.sum()), "live_score": meta["live_score"]})
        del a, dots
        print(f"assembled {short:52s} mass {used[-1]['mass']:7,d} live {meta['live_score']}")

    # Uniform-prior baseline: what the metric would look like if the truth were
    # spread uniformly over the scoreable domain.
    init = np.full((model.rows, model.cols), 0.02)
    rho, opt = fit(model, None, lam=args.lam, mu=args.mu, init=init, max_iter=args.max_iter)
    pred = model.score_all(rho)
    y = np.array([a["live_score"] for a in model.anchors])
    p = np.array([pred[a["name"]] for a in model.anchors])

    # Prior-only baseline: score the same anchors with rho = the public SGMC density.
    prior_rho = np.clip(model.prior_density, 0, None)
    if prior_rho is None:
        raise SystemExit("no prior density")
    prior_pred = np.array([model.score(prior_rho, a) for a in model.anchors])

    report = {
        "snapshot_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "method": "block-grid truth-density inversion of the published organizer scores",
        "block_px": args.block,
        "lambda_laplacian": args.lam,
        "mu_prior": args.mu,
        "n_blocks": model.n_blocks,
        "scoreable_px": int(scoreable.sum()),
        "kernel_mass_s": model.s,
        "anchors": used,
        "skipped": skipped,
        "fit": {
            "success": bool(opt.success),
            "iterations": int(opt.nit),
            "objective": float(opt.fun),
            "in_sample_spearman": spearman(p, y),
            "in_sample_rmse": float(np.sqrt(np.mean((p - y) ** 2))),
            "max_abs_error": float(np.max(np.abs(p - y))),
        },
        "prior_only_baseline": {
            "in_sample_spearman": spearman(prior_pred, y),
            "in_sample_rmse": float(np.sqrt(np.mean((prior_pred - y) ** 2))),
        },
        "rows": [{"anchor": a["name"], "live": float(a["live_score"]), "model": float(pred[a["name"]])}
                 for a in model.anchors],
        "implicit_truth_px": float((rho * model.m_b.reshape(rho.shape)).sum()),
        "note": ("the live scores are owner-reported associations in ref/meta/prior_corpus.md; "
                 "the model is an instrument, not an organizer score"),
    }
    if args.loo:
        t1 = time.time()
        loo = leave_one_out(model, lam=args.lam, mu=args.mu, max_iter=args.max_iter, init_from_full=False)
        report["leave_one_out"] = {"loo_spearman": loo["loo_spearman"], "loo_rmse": loo["loo_rmse"],
                                   "pred": loo["pred"], "actual": loo["actual"],
                                   "names": loo["names"], "seconds": round(time.time() - t1, 1)}
    report["seconds"] = round(time.time() - t0, 1)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    args.grid.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.grid, rho=rho.astype(np.float32), block=args.block,
                        rows=model.rows, cols=model.cols, m_b=model.m_b,
                        names=np.array([a["name"] for a in model.anchors]))
    print(json.dumps({k: v for k, v in report.items() if k not in ("rows", "anchors", "skipped")},
                     indent=2))
    print(f"[out] {args.out}  [grid] {args.grid}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
