#!/usr/bin/env python3
"""Fit and leave-one-out test instruments mapping measurable raster properties to the
owner-reported live scores of the anchor priors.

Motivation (measured, not assumed): this repository owns the *bytes* of the prior dot sets
whose public leaderboard scores are known, plus the provided catalogue and the public USGS
State Geologic Map Compilation (SGMC).  Everything the published metric needs for a
prediction raster can therefore be recomputed locally *except* the organizer's hidden
truth.  This script asks the question that decides whether local optimisation is
meaningful at all:

    can any locally measurable statistic rank the scored anchors the way the live
    leaderboard did, out of fold?

Four instruments are fitted on identical measurements and compared by leave-one-out
Spearman / RMSE:

``saturating``    S = TP/(0.2 TP + 0.2 M + 0.8 K),  TP = K(1 - exp(-a M w / K))
``powerlaw``      S = c * DTIproxy^p
``hidden``        S = lam*M*w / (0.2 lam*M*w + 0.2 M (1 - min(lam*w,1)) + 0.8 N)
                  (exact published metric algebra; lam = surrogate-to-hidden credit
                   ratio, N = hidden truth pixel count)
``block_linear``  S = c0 + sum_k c_k * (DTI of block k) over the 4x6 partition

A model is only usable as a promotion gate if its leave-one-out Spearman reaches the
pre-registered bar (default 0.80).  Nothing here is a competition score; ``live_score``
values are the owner-reported associations already recorded in ``ref/meta``.

    python scripts/fit_live_transfer.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import rasterio
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe40.measure import kernel_field, kernel_from_truth  # noqa: E402

DATA = ROOT / "data"
PRIOR = ROOT / "ref" / "prior"
META = ROOT / "ref" / "meta"
OUT = ROOT / "docs" / "data" / "live-transfer.json"
RADIUS_PX = 3.0
BLOCKS = (4, 6)
GUARD = 3


# --------------------------------------------------------------------------- #
# prediction functions (one per instrument) and their fitters
# --------------------------------------------------------------------------- #

def pred_saturating(params, mass, w):
    K, a = params
    tp = K * (1.0 - np.exp(-a * mass * w / K))
    return tp / (0.2 * tp + 0.2 * mass + 0.8 * K)


def fit_saturating(mass, w, y):
    res = optimize.minimize(lambda p: float(np.sum((pred_saturating(np.exp(p), mass, w) - y) ** 2)),
                            np.log([6_000.0, 1.0]), method="Nelder-Mead",
                            options=dict(xatol=1e-10, fatol=1e-14, maxiter=40_000))
    return np.exp(res.x)


def pred_powerlaw(params, proxy_dti):
    c, p = params
    return c * np.maximum(proxy_dti, 0.0) ** p


def fit_powerlaw(proxy_dti, y):
    sl, ic = np.polyfit(np.log(np.maximum(proxy_dti, 1e-9)), y, 1)
    return np.array([np.exp(ic), sl])


def pred_hidden(params, mass, w):
    lam, n_true = params
    t = lam * mass * w
    c = np.minimum(lam * w, 1.0) * mass
    den = 0.2 * t + 0.2 * (mass - c) + 0.8 * n_true
    return t / np.maximum(den, 1e-12)


def fit_hidden(mass, w, y):
    res = optimize.minimize(lambda p: float(np.sum((pred_hidden(np.exp(p), mass, w) - y) ** 2)),
                            np.log([1.0, 1e5]), method="Nelder-Mead",
                            options=dict(xatol=1e-10, fatol=1e-14, maxiter=60_000))
    return np.exp(res.x)


def pred_block(coef, per_block):
    return coef[0] + per_block @ coef[1:]


def fit_block(per_block, y):
    design = np.column_stack([np.ones(len(y)), per_block])
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    return coef


INSTRUMENTS = {
    "saturating": (fit_saturating, pred_saturating, ("mass", "w")),
    "powerlaw": (fit_powerlaw, pred_powerlaw, ("proxy_dti",)),
    "hidden": (fit_hidden, pred_hidden, ("mass", "w")),
    "block_linear": (fit_block, pred_block, ("per_block",)),
}


# --------------------------------------------------------------------------- #
# measurement
# --------------------------------------------------------------------------- #

def load_masks():
    with rasterio.open(DATA / "sample_submission.tif") as ds:
        footprint = np.isfinite(ds.read(1))
    with rasterio.open(DATA / "labels.tif") as ds:
        catalogue = ds.read(1) > 0
    with rasterio.open(DATA / "external" / "derived_sgmc_faults_100m_u8.tif") as ds:
        sgmc = ds.read(1) > 0
    return footprint, catalogue, sgmc


def anchor_scores() -> dict:
    """{id-or-filename -> owner-reported live score} from the two committed metadata files."""
    by_id: dict[str, float] = {}
    by_name: dict[str, float] = {}
    corpus = META / "prior_corpus.md"
    if corpus.exists():
        pattern = re.compile(
            r"\|\s*`([^`]+)`\s*\|\s*`([0-9a-f]+)`\s*\|[^|]*\|[^|]*\|[^|]*\|\s*([^|]+?)\s*\|")
        for line in corpus.read_text(encoding="utf-8").splitlines():
            m = pattern.match(line)
            if not m:
                continue
            raw = m.group(3).strip()
            if raw in ("", "—", "-"):
                continue
            try:
                by_id[m.group(2)] = float(raw)
            except ValueError:
                pass
    index = META / "prior_index.json"
    if index.exists():
        for short, entry in json.loads(index.read_text(encoding="utf-8")).items():
            if entry.get("live_score") is not None:
                by_name[Path(entry["path"]).name] = float(entry["live_score"])
                by_name[short] = float(entry["live_score"])
    return {"by_id": by_id, "by_name": by_name}


def sha256_16(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def block_rectangles(footprint, catalogue, truth):
    """Interior rectangles of the 4 x 6 partition, identical to measure.blocked_components.

    Because the guard band is >= the 3 px kernel radius, a dot inside a rectangle can only
    credit truth inside the same rectangle, so the rectangle DTI is a masked reduction of
    two global kernel fields (exact and fast).  Rectangles with fewer than 25 truth pixels
    are dropped.
    """
    h, w = footprint.shape
    y_edges = np.linspace(0, h, BLOCKS[0] + 1, dtype=int)
    x_edges = np.linspace(0, w, BLOCKS[1] + 1, dtype=int)
    rects = []
    for i in range(BLOCKS[0]):
        for j in range(BLOCKS[1]):
            rs = slice(int(y_edges[i]) + (GUARD if i > 0 else 0),
                       int(y_edges[i + 1]) - (GUARD if i < BLOCKS[0] - 1 else 0))
            cs = slice(int(x_edges[j]) + (GUARD if j > 0 else 0),
                       int(x_edges[j + 1]) - (GUARD if j < BLOCKS[1] - 1 else 0))
            t = truth[rs, cs] & footprint[rs, cs] & ~catalogue[rs, cs]
            if int(t.sum()) >= 25:
                rects.append((f"r{i}c{j}", rs, cs, int(t.sum())))
    return rects


def measure(path: Path, footprint, catalogue, truth, kernel_truth, rects) -> dict:
    with rasterio.open(path) as ds:
        a = ds.read(1)
    inside = np.isfinite(a) & (a > 0) & footprint
    scored = inside & ~catalogue
    n = int(scored.sum())
    out = {"file": path.name, "sha16": sha256_16(path), "mass": n,
           "on_catalogue": int((inside & catalogue).sum())}
    if n == 0:
        out.update({"w": 0.0, "proxy_tp": 0.0, "proxy_fp": 0.0, "proxy_dti": 0.0,
                    "per_block": [0.0] * len(rects)})
        return out
    out["w"] = float(kernel_truth[scored].mean())
    kf = kernel_field(scored.astype(np.float32))
    fp_weights = (1.0 - kernel_truth) * footprint
    tp = float(kf[truth & footprint].sum())
    fp = float((fp_weights * scored).sum())
    n_truth = int((truth & footprint).sum())
    out.update({"proxy_tp": tp, "proxy_fp": fp,
                "proxy_dti": float(tp / max(tp + 0.2 * fp + 0.8 * (n_truth - tp), 1e-9))})
    per_block = []
    for _, rs, cs, n_k in rects:
        tp_k = float(kf[rs, cs][truth[rs, cs] & footprint[rs, cs] & ~catalogue[rs, cs]].sum())
        fp_k = float((fp_weights[rs, cs] * scored[rs, cs]).sum())
        per_block.append(float(tp_k / max(tp_k + 0.2 * fp_k + 0.8 * (n_k - tp_k), 1e-9)))
    out["per_block"] = per_block
    return out


def spearman(a, b) -> float:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if a.size < 3 or np.all(a == a[0]) or np.all(b == b[0]):
        return float("nan")
    value = stats.spearmanr(a, b).statistic
    return float(value) if np.isfinite(value) else float("nan")


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--bar", type=float, default=0.80, help="pre-registered LOO Spearman bar")
    args = ap.parse_args()

    footprint, catalogue, sgmc = load_masks()
    truth = sgmc & ~catalogue & footprint
    kernel_truth = kernel_from_truth(truth, RADIUS_PX)
    rects = block_rectangles(footprint, catalogue, truth)
    scores = anchor_scores()
    print(f"footprint {int(footprint.sum()):,}  catalogue {int(catalogue.sum()):,}  "
          f"proxy truth {int(truth.sum()):,}  blocks {len(rects)}")

    rows = []
    for path in sorted(PRIOR.glob("*.tif")):
        m = measure(path, footprint, catalogue, truth, kernel_truth, rects)
        live = scores["by_name"].get(path.name) or scores["by_id"].get(m["sha16"])
        m["id"] = path.name if live is None else next(
            (k for k in scores["by_name"] if k.endswith(path.name)), path.name)
        m["live_score"] = live
        rows.append(m)
        print(f" {'*' if live is not None else ' '} {m['id'][:46]:46s} mass {m['mass']:8,d} "
              f"w {m['w']:.4f} proxyDTI {m['proxy_dti']:.4f} live {live}")

    scored = [r for r in rows if r["live_score"] is not None]
    print(f"\n{len(rows)} rasters measured, {len(scored)} with owner-reported live scores")
    y = np.array([r["live_score"] for r in scored], float)
    arrays = {"mass": np.array([r["mass"] for r in scored], float),
              "w": np.array([r["w"] for r in scored], float),
              "proxy_dti": np.array([r["proxy_dti"] for r in scored], float),
              "per_block": np.array([r["per_block"] for r in scored], float)}

    fits = {name: fit(*(arrays[k] for k in keys), y)
            for name, (fit, _, keys) in INSTRUMENTS.items()}
    report = {
        "n_rasters_measured": len(rows),
        "n_anchors_scored": len(scored),
        "proxy_truth": "USGS SGMC faults outside the provided catalogue, inside the sample footprint",
        "score_provenance": "owner-reported associations from ref/meta; not organizer receipts",
        "pre_registered_loo_bar": args.bar,
        "fits": {name: [float(v) for v in params] for name, params in fits.items()},
        "instruments": {},
    }

    for name, (fit, pred, keys) in INSTRUMENTS.items():
        in_sample = pred(fits[name], *(arrays[k] for k in keys))
        loo_pred = []
        for i in range(len(y)):
            keep = np.arange(len(y)) != i
            params = fit(*(arrays[k][keep] for k in keys), y[keep])
            loo_pred.append(float(pred(params, *(
                np.array([arrays[k][i]]) for k in keys))[0]))
        loo_pred = np.array(loo_pred, float)
        report["instruments"][name] = {
            "parameters": report["fits"][name],
            "in_sample_spearman": spearman(in_sample, y),
            "in_sample_rmse": float(np.sqrt(np.mean((in_sample - y) ** 2))),
            "loo_spearman": spearman(loo_pred, y),
            "loo_rmse": float(np.sqrt(np.mean((loo_pred - y) ** 2))),
            "loo_max_abs_error": float(np.abs(loo_pred - y).max()),
            "loo_pred": [float(v) for v in loo_pred],
            "passes_bar": bool(spearman(loo_pred, y) >= args.bar),
        }

    report["anchors"] = [{k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                          for k, v in r.items()} for r in scored]
    report["unscored_measured"] = [{"id": r["id"], "mass": r["mass"], "w": r["w"],
                                    "proxy_dti": r["proxy_dti"]} for r in rows
                                   if r["live_score"] is None]
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"\ninstrument comparison (leave-one-out over {len(scored)} scored anchors)")
    for name, blob in report["instruments"].items():
        print(f"  {name:13s} in-sample rho {blob['in_sample_spearman']:+.3f}  "
              f"LOO rho {blob['loo_spearman']:+.3f}  LOO RMSE {blob['loo_rmse']:.4f}  "
              f"max err {blob['loo_max_abs_error']:.4f}  "
              f"{'PASS' if blob['passes_bar'] else 'fail'}")
    print(f"[out] {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
