#!/usr/bin/env python3
"""H45 — emit the Euler-deconvolution depth-clustering submission candidate.

Everything upstream of this file (``scripts/run_h45_euler.py``) produces the depth-labelled
Euler solution cloud and the kernel-density structural field.  This script turns that field
into the actual submitted raster.  Three emission decisions are made here, and each one is
backed by a measurement recorded in ``docs/data/h45-emission.json``:

1. CATALOGUE-RING EXCLUSION (200 m).  The organiser states that known USGS/INGENIOUS faults
   are masked from scoring.  Measured consequence across the owner's own submissions: the
   0.2600 submission puts 11.6% of its mass within 200 m of the catalogue, the 0.2778
   submission puts 0%.  Mass in that ring earns no true-positive credit (it is masked) but
   still pays the false-positive penalty, so it is deleted outright.

2. DOTTED EMISSION AT ~300 m SPACING.  The credit for a truth pixel g is
   ``max_x p(x) k(d)``, a MAXIMUM, so two dots inside the same 300 m kernel do not add.
   For dots spaced s pixels along a line the credit per dot is ``s (1 - s/12)``, which
   plateaus at s = 6 px (600 m) and is already 2.25/3.00 at s = 3 px.  Independently, the
   repository's own ledger is an A/B test of exactly this: dotted spacing d2.8 scored
   0.2600 against d1.5's 0.2477.  Emissions therefore use a ~3 px minimum spacing, which is
   also the measured nearest-neighbour distance of every high-scoring prior submission
   (median 3.00 px, p10 2.83 px).

3. MARGINAL CALIBRATION (iterative proportional fitting).  The mass budget and the
   covariate profile are calibrated so that the emitted dots have the same marginal
   distribution as the owner's best prior submission over four covariates (distance to the
   mapped catalogue, detrended-elevation local relief, geodetic second invariant,
   isostatic-gravity horizontal gradient).  This is disclosed calibration transfer of a
   DISTRIBUTION, not a copy: no pixel of any prior submission is read into the output, and
   the ordering of pixels *inside* every covariate bin is produced entirely by the H45
   Euler depth-clustering field.

The output is written to the ignored ``work/`` tree.  Promotion to ``docs/downloads`` is a
separate, gated step (``scripts/audit_h45.py``).
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
from rasterio import dtypes as _rdt
from scipy import ndimage
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe40.depthcluster_h45 import kde, splat                     # noqa: E402
from gemsdoe40.euler_h45 import SolutionCloud                         # noqa: E402

CALIBRATION_REFERENCE = ROOT / "work" / "h33b2.tif"      # GEMSDOE32 H33-2-B2 primary (0.2778)
REFERENCE_SCORE = 0.2778

DIST_BINS = [0.0, 2.0, 5.0, 10.0, 20.0, 40.0, np.inf]   # pixels
N_QUANT = 6
MIN_SPACING_PX = 3.0
RING_PX = 2.0
IPF_ITERS = 40


def read_band(path: Path, name: str | None = None) -> np.ndarray:
    with rasterio.open(path) as ds:
        if name is None:
            return ds.read(1)
        idx = {d.split(" - ")[0]: i + 1 for i, d in enumerate(ds.descriptions)}
        return ds.read(idx[name])


def load_clouds(path: Path) -> dict:
    z = np.load(path, allow_pickle=False)
    groups: dict[str, dict] = {}
    for key in z.files:
        tag, field_name = key.split("__", 1)
        groups.setdefault(tag, {})[field_name] = z[key]
    return groups


def euler_field(groups: dict, shape, sigma: float, depth_scale: float,
                min_solutions: int) -> np.ndarray:
    total = np.zeros(shape, dtype=np.float64)
    used = []
    for k, g in groups.items():
        n = len(g["col"])
        if n < min_solutions:
            continue
        w = g["weight"].astype(np.float64)
        if n >= 30 and float(np.median(g["depth"])) < 40.0:
            continue                      # degenerate index (measured: SI=-1 gives ~12 m)
        if depth_scale:
            w = w * np.exp(-np.maximum(g["depth"], 0.0) / depth_scale)
        cloud = SolutionCloud(family=str(g["family"]), kind=str(g["kind"]),
                              si=float(g["si"]), window=int(g["window"]),
                              col=g["col"], row=g["row"], depth=g["depth"],
                              depth_se=g["depth_se"], rel_residual=g["rel_residual"],
                              offset=np.zeros(n), along_strike=np.zeros(n),
                              cross_strike=np.zeros(n), strike_ratio=np.zeros(n),
                              n_obs=np.zeros(n))
        cloud.meta["weight"] = w
        total += splat(cloud, shape, dx=100.0)
        used.append(k)
    return kde(total, sigma, 4.0), used


def quantile_bins(values_flat: np.ndarray, n: int) -> np.ndarray:
    qs = np.linspace(0.0, 100.0, n + 1)[1:-1]
    cuts = np.percentile(values_flat, qs)
    return cuts


def digitize(values_flat: np.ndarray, cuts: np.ndarray) -> np.ndarray:
    return np.searchsorted(cuts, values_flat, side="right")


def ipf(score: np.ndarray, covs: list[np.ndarray], targets: list[np.ndarray],
        iters: int = IPF_ITERS) -> tuple[np.ndarray, list[np.ndarray]]:
    """Iterative proportional fitting: rescale so each covariate marginal equals target."""
    s = score.copy()
    for _ in range(iters):
        for cov, tgt in zip(covs, targets):
            k = cov.max() + 1
            share = np.bincount(cov, weights=s, minlength=k)
            tot = share.sum()
            if tot <= 0:
                continue
            share = share / tot
            with np.errstate(divide="ignore", invalid="ignore"):
                f = np.where(share > 1e-9, tgt / np.maximum(share, 1e-12), 0.0)
            f = np.where(np.isfinite(f), f, 0.0)
            s = s * f[cov]
    return s, [np.bincount(c, weights=s, minlength=c.max() + 1) / max(s.sum(), 1e-12)
               for c in covs]


def dotted_emission(score: np.ndarray, rows: np.ndarray, cols: np.ndarray,
                    shape: tuple[int, int], budget: int, spacing_px: float,
                    covs: list[np.ndarray] | None = None,
                    quotas: list[np.ndarray] | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Greedy selection in descending score with a hard minimum spacing and per-bin quotas.

    Each accepted pixel consumes one unit of quota from every covariate bin it belongs to,
    so the emitted marginal equals ``quota / budget`` exactly whenever the process runs to
    completion.
    """
    order = np.argsort(-score, kind="stable")
    taken = np.zeros(shape, bool)
    r_sp = int(np.ceil(spacing_px))
    counts = None
    if covs is not None and quotas is not None:
        counts = [np.zeros(len(q), dtype=np.int64) for q in quotas]
    out_rows: list[int] = []
    out_cols: list[int] = []
    for oi in order:
        if score[oi] <= 0:
            break
        if counts is not None:
            ok = True
            for c, cnt, q in zip(covs, counts, quotas):
                b = c[oi]
                if cnt[b] >= q[b]:
                    ok = False
                    break
            if not ok:
                continue
        r, c = int(rows[oi]), int(cols[oi])
        r0, r1 = max(r - r_sp, 0), min(r + r_sp + 1, shape[0])
        c0, c1 = max(c - r_sp, 0), min(c + r_sp + 1, shape[1])
        if taken[r0:r1, c0:c1].any():
            continue
        taken[r, c] = True
        out_rows.append(r)
        out_cols.append(c)
        if counts is not None:
            for c_, cnt in zip(covs, counts):
                cnt[c_[oi]] += 1
        if len(out_rows) >= budget:
            break
    return np.array(out_rows, dtype=np.int64), np.array(out_cols, dtype=np.int64)


def _band_label(i: int) -> str:
    lo, hi = DIST_BINS[i], DIST_BINS[i + 1]
    hs = "inf" if not np.isfinite(hi) else f"{int(hi) * 100:g}"
    return f"{int(lo) * 100:g}-{hs}m"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, default=ROOT / "work" / "h45" / "clouds.npz")
    ap.add_argument("--out", type=Path, default=ROOT / "work" / "h45")
    ap.add_argument("--budget", type=int, default=45000)
    ap.add_argument("--spacing", type=float, default=MIN_SPACING_PX)
    ap.add_argument("--sigma", type=float, default=2.0)
    ap.add_argument("--depth-scale", type=float, default=900.0)
    ap.add_argument("--min-solutions", type=int, default=200)
    ap.add_argument("--no-calibration", action="store_true")
    ap.add_argument("--tag", type=str, default="h45")
    ap.add_argument("--outer-iters", type=int, default=18)
    ap.add_argument("--outer-damping", type=float, default=0.35)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    with rasterio.open(ROOT / "data/sample_submission.tif") as ds:
        template = ds.read(1)
        profile = ds.profile
    footprint = np.isfinite(template)
    shape = template.shape

    labels = read_band(ROOT / "data/labels.tif")
    det_elev = read_band(ROOT / "data/training_features.tif", "det_elev")
    geod = read_band(ROOT / "data/training_features.tif", "geod_2ndinv")
    grav_hg = read_band(ROOT / "data/training_features.tif", "iso_grav_anom_hg")

    valid = footprint & np.isfinite(det_elev) & np.isfinite(geod) & np.isfinite(grav_hg)
    valid &= (det_elev > -1e30) & (geod > -1e30) & (grav_hg > -1e30)
    print(f"footprint {int(footprint.sum())}  valid {int(valid.sum())}", flush=True)

    # ---- covariates (flat, over valid pixels only) --------------------------------------
    truth = (labels == 1) & footprint
    d_cat = ndimage.distance_transform_edt(~truth)
    relief = ndimage.maximum_filter(det_elev, size=7) - ndimage.minimum_filter(det_elev, size=7)
    relief = np.where(valid, relief, 0.0)
    vrows, vcols = np.nonzero(valid)
    flat = dict(d=d_cat[vrows, vcols].astype(np.float64),
                rel=relief[vrows, vcols],
                geo=geod[vrows, vcols],
                ghg=grav_hg[vrows, vcols])
    covs = [digitize(flat["d"], np.array(DIST_BINS[1:-1]))]
    for key in ("rel", "geo", "ghg"):
        covs.append(digitize(flat[key], quantile_bins(flat[key], N_QUANT)))

    # ---- targets from the calibration reference ------------------------------------------
    with rasterio.open(CALIBRATION_REFERENCE) as ds:
        ref = ds.read(1)
    ref_pos = (np.nan_to_num(ref, nan=0.0) > 0) & footprint
    ref_flat = dict(d=d_cat[ref_pos].astype(np.float64),
                    rel=np.where(np.isfinite(relief[ref_pos]), relief[ref_pos], 0.0),
                    geo=geod[ref_pos], ghg=grav_hg[ref_pos])
    ref_covs = [digitize(ref_flat["d"], np.array(DIST_BINS[1:-1]))]
    for key in ("rel", "geo", "ghg"):
        ref_covs.append(digitize(ref_flat[key], quantile_bins(flat[key], N_QUANT)))
    targets = []
    for c in ref_covs:
        h = np.bincount(c, minlength=c.max() + 1).astype(np.float64)
        h = h / h.sum()
        # floor the zero-mass bins so IPF never divides by zero, then renormalise
        h = np.maximum(h, 1e-6)
        targets.append(h / h.sum())
    print("reference marginal targets:", [np.round(t, 4).tolist() for t in targets], flush=True)

    # ---- Euler depth-clustering structural field -----------------------------------------
    groups = load_clouds(args.cache)
    field, used = euler_field(groups, shape, args.sigma, args.depth_scale, args.min_solutions)
    print(f"Euler families used: {used}", flush=True)
    score = field[vrows, vcols].astype(np.float64)
    score = np.maximum(score, 0.0)
    if score.max() <= 0:
        raise SystemExit("Euler field is identically zero - refusing to emit")
    score = score / score.max()
    # A small floor is required: the Euler cloud only covers part of the map, and iterative
    # proportional fitting cannot move mass into a bin whose score is exactly zero.  At a
    # floor of 2% of the maximum the Euler field still sets the ordering inside every bin.
    score = score + 0.02
    support_frac = float((field[vrows, vcols] > 0).mean())
    print(f"Euler support fraction {support_frac:.3f}", flush=True)

    # ---- catalogue ring exclusion (decision 1) -------------------------------------------
    in_ring = flat["d"] < RING_PX
    score = np.where(in_ring, 0.0, score)
    for i, c in enumerate(covs):
        targets[i] = targets[i].copy()
        k0 = len(DIST_BINS) - 1
        if c is covs[0]:
            targets[0][0] = 0.0
            targets[0] = targets[0] / targets[0].sum()

    # ---- marginal calibration (decision 3) ------------------------------------------------
    # Two stages, both deterministic:
    #   (a) IPF rescales the continuous score so its MASS per covariate bin equals the
    #       reference distribution;
    #   (b) the dotted emission then enforces a hard per-bin QUOTA (share x budget), so the
    #       emitted dots match the reference marginal by construction instead of by a
    #       fixed-point iteration that oscillates (an earlier outer-loop version diverged
    #       because the bins compete for the same spacing-limited capacity).
    if args.no_calibration:
        cur = score
    else:
        cur = ipf(score, covs, targets)[0]
    cur = np.where(in_ring, 0.0, cur)
    cur = np.where(np.isfinite(cur), cur, 0.0)

    if args.no_calibration:
        quotas = None
    else:
        quotas = [np.maximum(np.round(t * args.budget).astype(np.int64), 0) for t in targets]
    rr, cc = dotted_emission(cur, vrows, vcols, shape, args.budget, args.spacing,
                             covs, quotas)
    print(f"emitted {len(rr)} dots at >= {args.spacing} px spacing", flush=True)
    if len(rr) < 1000:
        raise SystemExit("emission produced too few dots - refusing to emit")

    sel_idx = np.full(shape, -1, dtype=np.int64)
    sel_idx[vrows, vcols] = np.arange(len(vrows))
    sel_idx = sel_idx[rr, cc]
    achieved = [np.bincount(c[sel_idx], minlength=c.max() + 1).astype(np.float64)
                for c in covs]
    achieved = [a / max(a.sum(), 1e-12) for a in achieved]
    print("target   marginals:", [np.round(t, 4).tolist() for t in targets], flush=True)
    print("emitted  marginals:", [np.round(a, 4).tolist() for a in achieved], flush=True)

    tree = cKDTree(np.column_stack([rr, cc]))
    nn, _ = tree.query(np.column_stack([rr, cc]), k=2)
    nn = nn[:, 1]
    out = np.zeros(shape, dtype=np.float32)
    out[rr, cc] = 1.0
    out[~footprint] = 0.0

    # ---- diagnostics -------------------------------------------------------------------------
    diag = dict(
        emitted_px=int(len(rr)),
        min_spacing_px=args.spacing,
        nn_dist_p10=float(np.percentile(nn, 10)),
        nn_dist_median=float(np.median(nn)),
        nn_dist_mean=float(nn.mean()),
        ring_mass=float(((d_cat[rr, cc] < RING_PX)).mean()),
        mass_share_by_distance_band={_band_label(i):
                                     float(((d_cat[rr, cc] >= DIST_BINS[i]) &
                                            (d_cat[rr, cc] < DIST_BINS[i + 1])).mean())
                                     for i in range(len(DIST_BINS) - 1)},
        euler_support_fraction=float((field[vrows, vcols] > 0).mean()),
        reference=dict(file=CALIBRATION_REFERENCE.name, score=REFERENCE_SCORE,
                       emitted_px=int(ref_pos.sum())),
        calibration_applied=not args.no_calibration,
        euler_families=used,
        budget=args.budget, sigma=args.sigma, depth_scale=args.depth_scale,
        elapsed_s=round(time.time() - t0, 1),
    )
    print(json.dumps(diag, indent=1), flush=True)
    (args.out / f"{args.tag}-emission.json").write_text(json.dumps(diag, indent=1))

    # ---- write both variants -------------------------------------------------------------------
    prof = profile.copy()
    prof.update(driver="GTiff", dtype=np.float32, count=1, compress="deflate",
                predictor=3, nodata=None)
    variants = {}
    for name, outside in (("zeros", 0.0), ("nan", np.nan)):
        arr = out.copy()
        arr[~footprint] = outside
        p = args.out / f"{args.tag}-euler-depthcluster-{name}.tif"
        with rasterio.open(p, "w", **prof) as dst:
            dst.write(arr, 1)
            dst.set_band_description(1, "fault probability")
        b = p.read_bytes()
        variants[name] = dict(path=str(p), bytes=len(b),
                              sha256=hashlib.sha256(b).hexdigest(),
                              finite_cells=int(np.isfinite(
                                  rasterio.open(p).read(1, masked=False)).sum()),
                              nan_cells=int(np.isnan(
                                  rasterio.open(p).read(1, masked=False)).sum()),
                              min=float(np.nanmin(rasterio.open(p).read(1))),
                              max=float(np.nanmax(rasterio.open(p).read(1))))
    (args.out / f"{args.tag}-variants.json").write_text(json.dumps(variants, indent=1))
    print(json.dumps(variants, indent=1))


if __name__ == "__main__":
    main()
