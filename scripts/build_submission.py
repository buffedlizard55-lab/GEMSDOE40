#!/usr/bin/env python
"""ARCHIVED ONLY — opt-in historical reproduction, not current submission advice.

Build the submission artifacts from the Euler depth-clustered solution cloud.

Three artifacts are written, all from the same physics, all float32 with values in [0, 1] inside
the 5,167,373-pixel footprint of the sample submission and NaN outside it, EPSG:32611, 100 m,
3730 x 3292, transform (100, 0, 243350, 0, -100, 4508550):

  A  ``gems40-euler-si0-depthcluster-crossfamily``  *the brief-mandated artifact*.
     Cross-family (magnetic AND gravity) Euler SI = 0 depth-cluster density, reduced to its
     positive-mass crests so that the emission is a 1-px-wide continuous ridge field rather than
     a 2-D blur; values are the normalised cluster density, so tight/shallow clusters keep the
     highest values.  This is the artifact whose spatial pattern is compared, pixel by pixel,
     with every prior submission in ``scripts/novelty_audit.py``.

  B  ``gems40-euler-si0-depthcluster-crossfamily-continuous``
     The same field *without* crest reduction: the literal "continuous raster via kernel density
     estimation, normalised to [0, 1]" of the brief.  Kept for completeness and audit; it emits
     ~4 M pixels of small mass, which the metric charges alpha for, so it is not the recommended
     upload.

  C  ``gems40-euler-augmented-<incumbent>``
     The incumbent site-best emission plus the Euler crest pixels that are new and more than
     200 m from the catalogue.  This artifact is *not* presented as a novel pattern (it inherits
     a previous submission's pixels by construction); it exists only because the measured
     marginal efficiency of the Euler crest pixels (``scripts/nms_refine_test.py``) exceeds the
     metric's break-even, which makes it the only variant here with a measured reason to expect a
     higher score than the current site best.

Every receipt records the 12 portal checks, the sha256, the emitted mass and the source-layer
provenance.  Nothing is uploaded by this script.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import distance_transform_edt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gems40.grid import footprint_from_sample, write_submission
from gems40.stack import (GRAV_LAYER, MAG_LAYERS, crest, crossfamily, load_fields,
                          multiscale_rankmin)

BANDS = {"rtp": 2, "tmi": 14, "mag_anom": 1, "iso_grav_anom": 13}
TARGET_BUDGET = 20_000
MULTISCALE_WINDOWS = (10, 15, 20, 30)   # H40-4 conjunction; see data/evidence/multiscale.json
MULTISCALE_BUDGET = 60_000              # added px; largest argmax of projected live gain
                                        # that stays inside the LM validity domain (<120 k px)
                                        # and is positive in 4/4 blocked folds     # crest pixels kept in artifact A (validated in scripts/nms_refine_test.py)
AUGMENT_BUDGET = None      # None = every off-catalogue Euler crest pixel not already in the
#                            incumbent.  Chosen from data/evidence/augmentation_budget_sweep.json:
#                            the LM-calibrated gain rises monotonically (0.2716 -> 0.3085) and stays
#                            4/4-folds-positive up to the full set (23,056 new px, 60,710 total,
#                            inside the instrument's stated validity domain of < 120 k px).


def main() -> int:
    from gemsdoe40.legacy_guard import require_legacy_opt_in
    require_legacy_opt_in()
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="work")
    ap.add_argument("--window", type=int, default=10,
                    help="Euler window scale (px) whose fields are shipped")
    ap.add_argument("--data", default="data")
    ap.add_argument("--prior", default="/tmp/data/prior")
    ap.add_argument("--outdir", default="docs/downloads")
    ap.add_argument("--stamp", default="20261005")
    ap.add_argument("--incumbent", default="top_02778_h33b2_zeros.tif")
    args = ap.parse_args()

    foot = footprint_from_sample(Path(args.data) / "example_submission.tif")
    fields = load_fields(args.work, window=args.window,
                         needed=tuple(MAG_LAYERS) + (GRAV_LAYER,))
    S = crossfamily(fields, foot)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    receipts = []

    # ---------------------------------------------------------------- A: brief-mandated artifact
    c = crest(np.nan_to_num(S, nan=0.0), foot)
    cv = np.where(c, np.nan_to_num(S, nan=0.0), -np.inf)
    vals = np.sort(cv[foot][np.isfinite(cv[foot])])[::-1]
    thr = vals[min(TARGET_BUDGET, vals.size) - 1] if vals.size else 0.0
    A = np.where(c & (np.nan_to_num(S, nan=0.0) >= thr), np.nan_to_num(S, nan=0.0), 0.0)
    A = np.where(foot, A, np.nan)
    digest = hashlib.sha256(np.nan_to_num(A, nan=-1.0).astype(np.float32).tobytes()).hexdigest()[:8]
    nameA = f"gems40-euler-si0-depthcluster-crossfamily-{args.stamp}-{digest}.tif"
    recA = write_submission(outdir / nameA, A, foot, nan_outside=True)
    recA["design"] = dict(method="Euler deconvolution SI=0 (Reid et al. 1990) on magnetic and "
                                 "gravity layers; depth-aware cluster weighting; KDE; 1-px crests",
                          layers=[f"{k} (band {b})" for k, b in BANDS.items()],
                          family_rule="element-wise min of the magnetic-family mean and the gravity KDE",
                          budget_px=int((A > 0).sum()), target_budget=TARGET_BUDGET)
    recA["ip_reminder"] = ("UNIQUE PATTERN: measured Jaccard overlap with the best prior submission "
                           "is ~0.003 (see data/evidence/novelty.json).")
    receipts.append(recA)
    print(f"[A] {nameA}: {recA['emitted_pixels']} px, sha256 {recA['sha256'][:16]}")

    # ---------------------------------------------------------------- B: continuous KDE variant
    B = np.where(foot, S, np.nan)
    digest_b = hashlib.sha256(np.nan_to_num(B, nan=-1.0).astype(np.float32).tobytes()).hexdigest()[:8]
    nameB = f"gems40-euler-si0-depthcluster-crossfamily-continuous-{args.stamp}-{digest_b}.tif"
    recB = write_submission(outdir / nameB, B, foot, nan_outside=True)
    recB["design"] = dict(method="literal brief reading: continuous normalised KDE, no crest reduction",
                          budget_px=int((B > 0).sum()))
    recB["warning"] = ("~4 M pixels of small mass; the metric charges alpha per emitted pixel, so "
                       "this variant scores far below artifact A. Kept for audit only.")
    receipts.append(recB)
    print(f"[B] {nameB}: {recB['emitted_pixels']} px, sha256 {recB['sha256'][:16]}")

    # ---------------------------------------------------------------- C: validated augmentation
    inc_path = Path(args.prior) / args.incumbent
    if inc_path.exists():
        import rasterio
        with rasterio.open(inc_path) as ds:
            inc = (np.isfinite(ds.read(1)) & (ds.read(1) > 0))
        d_cat = d_cat_ms = distance_transform_edt(~_read_labels(args.data))
        if AUGMENT_BUDGET is None:
            add = c & (d_cat > 2) & ~inc
        else:
            thr_c = vals[min(AUGMENT_BUDGET, vals.size) - 1] if vals.size else 0.0
            add = c & (np.nan_to_num(S, nan=0.0) >= thr_c) & ~inc & (d_cat > 2)
        C = np.where(foot, np.maximum(np.nan_to_num(A, nan=0.0), add.astype(float)), np.nan)
        C = np.where(foot & (inc | add), np.where(inc, 1.0, np.nan_to_num(S)), np.nan)
        C = np.where(foot, np.nan_to_num(C, nan=0.0), np.nan)
        digest_c = hashlib.sha256(np.nan_to_num(C, nan=-1.0).astype(np.float32).tobytes()).hexdigest()[:8]
        nameC = f"gems40-euler-augmented-incumbent-{args.stamp}-{digest_c}.tif"
        recC = write_submission(outdir / nameC, C, foot, nan_outside=True)
        recC["design"] = dict(method=f"prior artifact '{args.incumbent}' union Euler crest pixels "
                                     f"(top {AUGMENT_BUDGET} of the cross-family field, >200 m from "
                                     "the catalogue)",
                              added_px=int(add.sum()),
                              jaccard_with_incumbent=float((C > 0)[inc].mean()))
        recC["not_unique_warning"] = ("contains a previous submission's pixels; documented as an "
                                      "augmentation, not as a novel pattern")
        recC["validation"] = json.loads((Path("data/evidence/nms_refine.json")).read_text())["rows"] \
            if Path("data/evidence/nms_refine.json").exists() else None
        receipts.append(recC)
        print(f"[C] {nameC}: {recC['emitted_pixels']} px (+{int(add.sum())} new), sha256 {recC['sha256'][:16]}")

    # ------------------------------------------------- D: multi-scale (H40-4) augmentation
    # Measured in scripts/multiscale_test.py (data/evidence/multiscale.json): a cross-family
    # rank-min conjunction over four Euler window scales, crest-reduced and capped at the budget,
    # emits 97,654 px in total (inside the instrument's <120k validity domain), earns
    # 1.18 x the metric's break-even per added pixel, and beats the incumbent in 4/4 blocked folds
    # with a projected live gain of +0.029 -- 2.6x the single-scale rule it supersedes.
    if inc_path.exists():
        import rasterio as _rio
        ms_fields = {w: crossfamily(load_fields(args.work, window=w,
                                                needed=tuple(MAG_LAYERS) + (GRAV_LAYER,)), foot)
                     for w in MULTISCALE_WINDOWS}
        R = multiscale_rankmin(ms_fields, foot)
        with _rio.open(inc_path) as ds:
            inc_ms = (np.isfinite(ds.read(1)) & (ds.read(1) > 0))
        pool = crest(np.nan_to_num(R, nan=0.0), foot) & ~inc_ms & (d_cat_ms > 2)
        flat = np.where(pool, np.nan_to_num(R, nan=0.0), -1.0)
        k = min(MULTISCALE_BUDGET, int(pool.sum()))
        thr = np.partition(flat[foot], int(foot.sum()) - k)[int(foot.sum()) - k]
        add_ms = pool & (np.nan_to_num(R, nan=0.0) >= thr)
        D = np.where(foot & (inc_ms | add_ms), np.where(inc_ms, 1.0, np.nan_to_num(R)), np.nan)
        D = np.where(foot, np.nan_to_num(D, nan=0.0), np.nan)
        digest_d = hashlib.sha256(np.nan_to_num(D, nan=-1.0).astype(np.float32).tobytes()).hexdigest()[:8]
        nameD = f"gems40-euler-multiscale-augmented-incumbent-{args.stamp}-{digest_d}.tif"
        recD = write_submission(outdir / nameD, D, foot, nan_outside=True)
        recD["design"] = dict(method=f"prior artifact '{args.incumbent}' union the multi-scale "
                                     f"rank-min Euler conjunction crests, top {k} new pixels, "
                                     f">200 m from the catalogue",
                              windows=list(MULTISCALE_WINDOWS), added_px=int(add_ms.sum()),
                              evidence="data/evidence/multiscale.json")
        recD["not_unique_warning"] = ("contains a previous submission's pixels; documented as an "
                                      "augmentation, not as a novel pattern")
        receipts.append(recD)
        print(f"[D] {nameD}: {recD['emitted_pixels']} px (+{int(add_ms.sum())} new), "
              f"sha256 {recD['sha256'][:16]}")

    # validator-safety companions: identical bytes inside the footprint, 0 (not NaN) outside.
    # The sample submission itself is NaN outside, but DrivenData's own validator rejects any
    # value it considers out of range, and a submission that failed with
    # "Predicted values must be in range [0, 1]" is exactly the symptom of NaN inside the
    # footprint -- these companions remove every possible ambiguity.
    for rec in receipts:
        if "continuous" in rec["filename"]:
            continue          # audit-only variant; the 13 MB zeros companion is not needed
        src = outdir / rec["filename"]
        with rasterio.open(src) as ds:
            arr = ds.read(1)
        zeros_name = rec["filename"].replace(".tif", "-zeros.tif")
        with rasterio.open(src) as ds:
            prof = ds.profile
        prof.update(nodata=None)
        z = np.where(np.isnan(arr), 0.0, arr).astype("float32")
        with rasterio.open(outdir / zeros_name, "w", **prof) as ds:
            ds.write(z, 1)
        rec["zeros_companion"] = zeros_name
        rec["zeros_companion_sha256"] = hashlib.sha256((outdir / zeros_name).read_bytes()).hexdigest()
        with rasterio.open(outdir / zeros_name) as ds:
            zz = ds.read(1)
        rec["zeros_companion_range_ok"] = bool(np.isfinite(zz).all() and (zz >= 0).all() and (zz <= 1).all())

    manifest = outdir / "manifest.json"
    manifest.write_text(json.dumps(receipts, indent=2))
    print(f"[out] {manifest}")
    return 0


def _read_labels(data_dir: str) -> np.ndarray:
    """Catalogue fault mask (USGS Quaternary + INGENIOUS), the population the scorer masks out."""
    import rasterio
    with rasterio.open(Path(data_dir) / "existing_faults.tif") as ds:
        return (ds.read(1) > 0)


if __name__ == "__main__":
    raise SystemExit(main())
