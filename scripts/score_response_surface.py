#!/usr/bin/env python3
"""Measure every cached prior prediction raster and fit the owner-reported score response.

Purpose
-------
The organiser's hidden test labels are unavailable to us.  What *is* available is a set of
owner-reported public-leaderboard scores for a subset of the public sibling-repository
predictions, plus the rasters themselves (pinned by git blob SHA-1 in
``docs/data/prior-inventory-20261006.json`` and cached under ``data/prior/``).

This script turns that into a measured response surface:

  1. For every cached raster that can be matched to an owner-reported score, compute
     purely geometric descriptors that do **not** require the hidden labels
     (emitted mass, positive-pixel count, hit counts against the public USGS/INGENIOUS
     label raster at several radii, hit counts against the independent SGMC-derived
     catalogue, spatial dispersion, value-distribution moments).
  2. Report the resulting table and the rank correlations between each descriptor and
     the reported score.

This is an *observational* regression over a small, family-correlated sample.  It is
evidence about which emission statistics accompany higher reported scores; it is not a
causal model, not an organiser receipt, and not a promise about any future score.
No hidden label, no leaderboard credential and no network call is used.

Usage
-----
    . .venv/bin/activate
    python scripts/score_response_surface.py [--cache data/prior] [--out work/]
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import re

import numpy as np
import rasterio
from scipy import ndimage
from scipy.ndimage import distance_transform_edt
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]

# --------------------------------------------------------------------------------------
# Owner-reported scores.  Two independent transcriptions of the same public statements are
# merged: the repository score ledger (ref/gems32/score-ledger.csv) and the retained user
# brief.  Where both exist they agree; the brief adds the GEMSDOE27-33 rows.
# provenance: "user-reported; not independently authenticated" for every entry.
# --------------------------------------------------------------------------------------
BRIEF_SCORES: dict[str, float] = {
    # GEMSDOE28
    "h27-4-r1-solo-d2-8-20261003-8acb75e1f2cc-nan": 0.2708,
    "h32-1-prethin-tip-euler-d2-8-20261003-31e35eee884e-nan": 0.2649,
    # GEMSDOE31
    "h27-4-solo-d28-20261004-8acb75e1-nan": 0.2708,
    # GEMSDOE32
    "h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros": 0.2778,
    # GEMSDOE33
    "h33d-analog-tip-stepover-r30-20261004-cb490425926e": 0.2632,
}

SUFFIX_RE = re.compile(r"[-_.](nan|zeros|allfinite|nanoutside|nan-outside|candidate|zero)$", re.I)


def normalise(name: str) -> str:
    """Strip the .tif extension and the outside-footprint encoding suffix."""
    name = re.sub(r"\.tif$", "", name, flags=re.I)
    prev = None
    while prev != name:
        prev = name
        name = SUFFIX_RE.sub("", name)
    return name


def load_scores(ledger: Path) -> dict[str, float]:
    scores: dict[str, float] = {}
    for row in csv.DictReader(ledger.open()):
        site = (row.get("site") or "").lower()
        if "drivendata" in site:          # leaderboard rows, not owner submissions
            continue
        try:
            value = float(row["score"])
        except (TypeError, ValueError):
            continue
        scores[row["submission"].strip()] = value
    for name, value in BRIEF_SCORES.items():       # brief adds rows absent from the ledger
        scores.setdefault(name, value)
    return scores


def match_score(filename: str, scores: dict[str, float]) -> tuple[float, str, str] | None:
    """Match a cached raster to an owner-reported score key.

    Only *token-boundary* containment is accepted.  Generic file names such as
    ``submission.tif`` must not match ``gems-submission-2026...`` simply because the
    word "submission" appears in both; requiring a ``-`` delimited boundary and a
    minimum key length keeps the response surface from being polluted by
    mis-attributed scores.  A leading repository tag in the file name
    (``gemsdoe23-``, ``gems10-`` ...) is tolerated, since sibling repos rename the
    same artifact with a repo prefix.
    """
    stem = normalise(Path(filename).name)
    best: tuple[float, str, str] | None = None
    for key, value in scores.items():
        k = normalise(key)
        if len(k) < 12:
            continue
        kind = None
        if k == stem:
            kind = "exact"
        elif stem.endswith("-" + k):          # repo-prefixed copy of the same artifact
            kind = "repo-prefix"
        elif k.endswith("-" + stem):
            kind = "key-suffix"
        elif ("-" + k + "-") in ("-" + stem + "-"):
            kind = "token-substring"
        if kind is None:
            continue
        if best is None or len(k) > len(normalise(best[1])):
            best = (value, key, kind)
    return best


# --------------------------------------------------------------------------------------
# Reference geometry (all hash-pinned public inputs; never the hidden labels)
# --------------------------------------------------------------------------------------
def load_references():
    with rasterio.open(ROOT / "data/sample_submission.tif") as ds:
        template = ds.read(1)
    footprint = np.isfinite(template)
    with rasterio.open(ROOT / "data/labels.tif") as ds:
        raw = ds.read(1)
        nodata = ds.nodata
    known = (raw == 1) & footprint if nodata is not None else (raw == 1) & footprint
    with rasterio.open(ROOT / "data/external/derived_sgmc_faults_100m_u8.tif") as ds:
        sg = ds.read(1)
    sgmc = (sg == 1) & footprint
    return footprint, known, sgmc


def tri_credit(mask: np.ndarray, pred: np.ndarray, radius: float = 3.0) -> tuple[float, float, int]:
    """Official distance-weighted TP/FP terms for a soft prediction against ``mask``.

    Mirrors src/gems40/metric.py: TP_w uses a per-truth-pixel max, FP_w sums
    p * (1 - max_g k(d)) over every predicted pixel.
    """
    offsets = [(dy, dx, max(1.0 - math.hypot(dy, dx) / radius, 0.0))
               for dy in range(-3, 4) for dx in range(-3, 4)
               if max(1.0 - math.hypot(dy, dx) / radius, 0.0) > 0.0]
    yy, xx = np.nonzero(mask)
    credit = np.zeros(yy.size, dtype=np.float64)
    if yy.size:
        for dy, dx, k in offsets:
            ny, nx = yy + dy, xx + dx
            ok = (ny >= 0) & (ny < pred.shape[0]) & (nx >= 0) & (nx < pred.shape[1])
            credit[ok] = np.maximum(credit[ok], pred[ny[ok], nx[ok]] * k)
    tp = float(credit.sum())
    near = np.maximum(1.0 - distance_transform_edt(~mask) / radius, 0.0) if mask.any() else np.zeros_like(pred)
    fp = float((pred * (1.0 - near)).sum())
    return tp, fp, int(yy.size)


def describe(path: Path, footprint, known, sgmc, d_known, d_sgmc) -> dict | None:
    try:
        with rasterio.open(path) as ds:
            arr = ds.read(1)
    except Exception:
        return None
    if arr.shape != footprint.shape:
        return None
    a = np.asarray(arr, dtype=np.float64)
    inside = footprint & np.isfinite(a)
    if inside.sum() == 0:
        return None
    p = np.where(inside, a, 0.0)
    p = np.where(p > 0, p, 0.0)
    mass = float(p.sum())
    if mass <= 0:
        return None
    pos = int((p > 0).sum())
    out = {
        "grid_matches": True,
        "finite_inside": int((footprint & np.isfinite(a)).sum()),
        "nan_inside": int((footprint & ~np.isfinite(a)).sum()),
        "positive_px": pos,
        "mass": mass,
        "max": float(p.max()),
        "mean_positive": float(p[p > 0].mean()),
        "distinct_values": int(np.unique(np.round(p[p > 0], 6)).size),
    }
    # mass-distribution shape
    q = np.quantile(p[p > 0], [0.1, 0.5, 0.9]) if pos else np.zeros(3)
    out["q10"], out["q50"], out["q90"] = (float(v) for v in q)
    out["mass_fraction_top10pct_positives"] = float(
        np.sort(p[p > 0])[int(0.9 * pos):].sum() / mass) if pos else 0.0
    # hits against the two public catalogues
    for radius, tag in ((1, "100m"), (2, "200m"), (3, "300m")):
        out[f"pos_within_{tag}_known"] = int(((d_known <= radius) & (p > 0)).sum())
    out["pos_within_300m_sgmc"] = int(((d_sgmc <= 3) & (p > 0)).sum())
    # mass inside the 300 m catalogue halo (the region the organiser masks / the owner prunes)
    halo = d_known <= 3
    out["mass_fraction_in_known_halo300m"] = float(p[halo].sum() / mass)
    out["mass_fraction_beyond_1km_known"] = float(p[d_known > 10].sum() / mass)
    # dispersion: median nearest-neighbour spacing of the emitted support (px).
    # A KD-tree keeps this O(n log n) in time and O(n) in memory; the earlier
    # blockwise brute force allocated a (block x n) distance matrix and exhausted RAM.
    if 1 < pos <= 600000:
        from scipy.spatial import cKDTree
        pts = np.column_stack(np.nonzero(p > 0)).astype(np.float32)
        step = max(1, pos // 40000)
        sample = pts[::step]
        distances, _ = cKDTree(pts).query(sample, k=2, workers=1)
        # k=2 returns self (0) first; the second column is the nearest *other* support pixel.
        out["mean_nn_spacing_px"] = float(np.median(distances[:, 1]))
    else:
        out["mean_nn_spacing_px"] = float("nan")
    # official DTI terms against each public catalogue (descriptive only)
    for tag, mask in (("known", known), ("sgmc", sgmc)):
        tp, fp, n_g = tri_credit(mask, p)
        out[f"tp_w_{tag}"] = tp
        out[f"fp_w_{tag}"] = fp
        out[f"credit_per_unit_mass_{tag}"] = tp / mass
        out[f"dti_proxy_{tag}"] = tp / (tp + 0.2 * fp + 0.8 * n_g + 1e-12)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, default=ROOT / "data/prior")
    parser.add_argument("--inventory", type=Path, default=ROOT / "docs/data/prior-inventory-20261006.json")
    parser.add_argument("--out", type=Path, default=ROOT / "work")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    scores = load_scores(ROOT / "ref/gems32/score-ledger.csv")
    inventory = json.loads(args.inventory.read_text())
    footprint, known, sgmc = load_references()
    d_known = distance_transform_edt(~known)
    d_sgmc = distance_transform_edt(~sgmc)

    rows = []
    unmatched = []
    for entry in inventory["unique_rasters"]:
        blob = entry["git_blob_sha"]
        cached = args.cache / f"{blob}.tif"
        if not cached.exists():
            continue
        hit = None
        for artifact in entry["artifacts"]:
            hit = match_score(Path(artifact["path"]).name, scores)
            if hit:
                break
        if not hit:
            unmatched.append(entry["artifacts"][0]["path"])
            continue
        value, key, kind = hit
        desc = describe(cached, footprint, known, sgmc, d_known, d_sgmc)
        if desc is None:
            continue
        rows.append(dict(blob=blob, repo=entry["artifacts"][0]["repo"],
                         path=entry["artifacts"][0]["path"], score=value,
                         score_key=key, match=kind, **desc))

    # One row per distinct submission *content*: the -nan / -zeros / -allfinite twins carry
    # identical pixel values and only differ in how the outside-footprint is encoded, so
    # keeping both would double-weight a single experiment in the correlation.
    seen: dict[tuple, dict] = {}
    for r in rows:
        fingerprint = (r["score_key"], r["positive_px"], round(r["mass"], 3))
        seen.setdefault(fingerprint, r)
    rows = list(seen.values())
    rows.sort(key=lambda r: -r["score"])
    payload = {
        "evidence_class": "observational fit over owner-reported public scores; not an organiser receipt",
        "n_scored_rasters": len(rows),
        "n_unmatched_cached": len(unmatched),
        "score_provenance": "user-reported public leaderboard values; not independently authenticated",
        "rows": rows,
        "correlations": {},
    }
    keys = ["mass", "positive_px", "mean_positive", "distinct_values",
            "pos_within_100m_known", "pos_within_200m_known", "pos_within_300m_known",
            "pos_within_300m_sgmc", "mass_fraction_in_known_halo300m",
            "mass_fraction_beyond_1km_known", "mean_nn_spacing_px",
            "credit_per_unit_mass_known", "credit_per_unit_mass_sgmc",
            "dti_proxy_known", "dti_proxy_sgmc"]
    vals = {k: np.array([r[k] for r in rows], dtype=float) for k in keys}
    y = np.array([r["score"] for r in rows], dtype=float)
    for k in keys:
        good = np.isfinite(vals[k])
        if good.sum() >= 5 and np.unique(vals[k][good]).size >= 3:
            rho, p = spearmanr(vals[k][good], y[good])
            payload["correlations"][k] = {"spearman_rho": float(rho), "p_value": float(p),
                                          "n": int(good.sum())}
    (args.out / "score_response_surface.json").write_text(json.dumps(payload, indent=1))

    print(f"matched {len(rows)} cached rasters to owner-reported scores "
          f"({len(unmatched)} cached rasters unmatched)")
    header = f"{'score':>7} {'repo':11s} {'px':>7} {'mass':>9} {'<=200m':>7} {'<=300m':>7} {'sgmc300':>7} {'cpu_kn':>7} {'cpu_sg':>7} {'nnsp':>5} file"
    print(header)
    for r in rows:
        print(f"{r['score']:7.4f} {r['repo']:11s} {r['positive_px']:7d} {r['mass']:9.0f} "
              f"{r['pos_within_200m_known']:7d} {r['pos_within_300m_known']:7d} "
              f"{r['pos_within_300m_sgmc']:7d} {r['credit_per_unit_mass_known']:7.4f} "
              f"{r['credit_per_unit_mass_sgmc']:7.4f} {r['mean_nn_spacing_px']:5.1f} "
              f"{Path(r['path']).name}")
    print()
    print("Spearman correlation with the reported score:")
    for k, v in sorted(payload["correlations"].items(), key=lambda kv: -abs(kv[1]["spearman_rho"])):
        print(f"  {k:34s} rho={v['spearman_rho']:+.3f}  p={v['p_value']:.4g}  n={v['n']}")


if __name__ == "__main__":
    main()
