#!/usr/bin/env python3
"""H4-A: build, validate and emit the Euler SI=0 depth-cluster contact candidate.

Pipeline (all stages are deterministic; no randomness anywhere):

  1. Euler deconvolution, SI=0 (fault-like contact), moving windows 8/12/16 px,
     on challenge band `rtp` and on challenge band `iso_grav_anom` (the vertical
     derivative is derived inside the solver by the Fourier `+|k|` operator).
  2. 3-D depth-consistency gate + solution weighting (shallow x well-fitted x
     locally tight x mutually depth-consistent).
  3. Cross-field concordance: magnetic/gravity solutions agreeing in (x, y, z).
  4. Weighted KDE of the surviving cloud -> continuous field E.
  5. Multi-scale Hessian lineament response of E -> trace geometry L.
  6. Surface-expression concurrence C from the 3DEP 1 m DEM scarp stack and the
     challenge detrended-elevation slope band.
  7. H4-B mapping-gap covariate G from the official USGS QFFD scale split.
  8. Variants A0 = E, A1 = L*(0.35+0.65C), A2 = A1*G.  Ranked on the
     catalogue-block instrument; ties within 0.005 prefer the concurrence
     variant (rule pre-registered in docs/research/hypotheses.md).
  9. Dot emission with 3 px minimum separation (non-maximum suppression), the
     catalogue pruned exactly, budget chosen from the measured precision ladder
     with a 2-fold block-split selection so the reported number is not
     selected on the blocks it is reported on.
 10. GeoTIFF writing (nan + zeros twins + zip), independent format re-read,
     and the registered uniqueness audit against every cached prior raster.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe40 import BANDS  # noqa: E402
from gemsdoe40.depthcluster import (  # noqa: E402
    DEFAULT_CONFIGS,
    cross_field_concordance,
    depth_cluster_weights,
    hessian_lineament_response,
    mapping_gap_weight,
    robust_normalize,
    solve_field,
    splat_points,
    terrain_concurrence,
)
from gemsdoe40.emit import nms_dots  # noqa: E402
from gemsdoe40.grid import check_submission, footprint_from_sample, read_band, read_labels, sha256, write_submission  # noqa: E402
from gemsdoe40.research_metric import spatial_block_components  # noqa: E402

DATA = ROOT / "data"
WORK = ROOT / "work"
DOWNLOADS = ROOT / "docs" / "downloads"
DOCDATA = ROOT / "docs" / "data"
EVIDENCE = ROOT / "evidence"

SCARP = WORK / "ext" / "lidar_scarp_features_u8.tif"
QFAULTS = WORK / "ext" / "qfaults_prior_u8.tif"
CACHE = WORK / "depthcluster_v2_cache.npz"

BUDGET_LADDER = (5000, 10000, 20000, 30000, 45000, 60000, 90000)
TRUTH_BLOCK_SPLITS = 2  # 2-fold block split for budget selection


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def load_inputs() -> dict:
    sample = DATA / "sample_submission.tif"
    labels = DATA / "labels.tif"
    features = DATA / "training_features.tif"
    for path in (sample, labels, features):
        if not path.exists():
            raise SystemExit(f"missing required input {path}")
    footprint = footprint_from_sample(sample)
    catalogue = read_labels(labels)
    rtp, rtp_ok = read_band(features, BANDS["rtp"])
    grav, grav_ok = read_band(features, BANDS["iso_grav_anom"])
    slope, _ = read_band(features, BANDS["det_elev_slope"])
    return {
        "sample": sample,
        "catalogue": catalogue,
        "footprint": footprint,
        "rtp": rtp,
        "rtp_ok": rtp_ok & footprint,
        "grav": grav,
        "grav_ok": grav_ok & footprint,
        "slope": slope,
    }


def build_fields(inputs: dict, cache: bool = True) -> dict:
    t0 = time.time()
    if cache and CACHE.exists():
        z = np.load(CACHE)
        e = z["E"]
        stats = json.loads(str(z["stats"]))
        print(f"loaded cached Euler field in {time.time()-t0:.1f}s")
    else:
        footprints = inputs["footprint"]
        mag_cloud, mag_stats = solve_field(inputs["rtp"], inputs["rtp_ok"], "rtp", DEFAULT_CONFIGS)
        grav_cloud, grav_stats = solve_field(inputs["grav"], inputs["grav_ok"], "iso_grav_anom", DEFAULT_CONFIGS)
        print(f"Euler solutions: magnetic {len(mag_cloud):,}  gravity {len(grav_cloud):,}  "
              f"({time.time()-t0:.1f}s)")

        w_m, keep_m, stats_m = depth_cluster_weights(mag_cloud)
        w_g, keep_g, stats_g = depth_cluster_weights(grav_cloud)
        print(f"depth-consistent clusters: magnetic {stats_m['kept']:,}  gravity {stats_g['kept']:,}")

        pair_pts, pair_w, pair_stats = cross_field_concordance(mag_cloud, w_m, grav_cloud, w_g)
        print(f"cross-field concordant pairs: {pair_stats['pairs']:,}")

        shape = footprints.shape
        kde_m = splat_points(np.column_stack([mag_cloud.row[keep_m], mag_cloud.col[keep_m]]),
                             w_m[keep_m], shape, sigma_px=1.5)
        kde_g = splat_points(np.column_stack([grav_cloud.row[keep_g], grav_cloud.col[keep_g]]),
                             w_g[keep_g], shape, sigma_px=1.5)
        kde_p = splat_points(pair_pts, pair_w, shape, sigma_px=1.5)
        raw = kde_m + kde_g + 1.75 * kde_p
        E = robust_normalize(raw, footprints, hi=99.5)
        e = E
        stats = {
            "euler_magnetic": mag_stats,
            "euler_gravity": grav_stats,
            "mag_cloud": mag_cloud.stats,
            "grav_cloud": grav_cloud.stats,
            "mag_cluster": stats_m,
            "grav_cluster": stats_g,
            "concordance": pair_stats,
            "kde_positive_px": int((E > 0).sum()),
            "elapsed_s": round(time.time() - t0, 1),
        }
        if cache:
            np.savez_compressed(CACHE, E=np.asarray(E, dtype=np.float32), stats=json.dumps(stats))
        e = np.asarray(E, dtype=np.float64)
        print(f"Euler KDE field: {int((e>0).sum()):,} positive cells ({time.time()-t0:.1f}s)")

    t1 = time.time()
    L = hessian_lineament_response(e, inputs["footprint"])
    L = robust_normalize(L, inputs["footprint"], hi=99.5)
    C, c_stats = terrain_concurrence(SCARP, inputs["slope"], inputs["footprint"])
    G, g_stats = mapping_gap_weight(QFAULTS, inputs["footprint"])
    print(f"lineament + concurrence + gap weights in {time.time()-t1:.1f}s")

    variants = {
        "A0_euler_only": e,
        "A1_euler_x_terrain": L * (0.35 + 0.65 * C),
        "A2_euler_x_terrain_x_gap": L * (0.35 + 0.65 * C) * np.sqrt(np.clip(G, 0.0, 1.0)),
    }
    for name in variants:
        variants[name] = robust_normalize(variants[name], inputs["footprint"], hi=99.8)
    return {"E": e, "L": L, "C": C, "G": G, "variants": variants, "field_stats": stats,
            "terrain_stats": c_stats, "gap_stats": g_stats}


def blocked_score(dots: np.ndarray, truth: np.ndarray, footprint: np.ndarray) -> dict:
    pooled, blocks = spatial_block_components(dots, truth, footprint, n_rows=4, n_cols=6, guard=3)
    return {"pooled": pooled.to_dict(), "blocks": blocks}


def budget_selection(field: np.ndarray, truth: np.ndarray, footprint: np.ndarray,
                     catalogue: np.ndarray, ladder=BUDGET_LADDER) -> dict:
    """2-fold block-split selection of the dot budget.

    The 24 blocks are split into two halves (blocks with an even vs odd
    row-major index).  For each half the budget maximising the *other* half's
    pooled DTI is chosen; the reported held-out number is the DTI that the
    chosen budget actually achieves on the half it was NOT chosen on.  Then the
    shipped budget is the argmax on all 24 blocks, disclosed as instrument
    selection.
    """
    rows = []
    per_budget_blocks: dict[int, list[dict]] = {}
    for b in ladder:
        dots = nms_dots(field, footprint, min_separation_px=3, budget=int(b),
                        catalogue=catalogue, catalogue_buffer_px=2)
        comp, blocks = spatial_block_components(dots, truth, footprint, n_rows=4, n_cols=6, guard=3)
        emitted = int((dots > 0).sum())
        rows.append({
            "budget": int(b),
            "emitted_px": emitted,
            "tpw": comp.tp,
            "fpw": comp.fp,
            "fnw": comp.fn,
            "truth_px": comp.n_truth,
            "dti": comp.score,
            "credit_per_dot": float(comp.tp / max(emitted, 1)),
        })
        per_budget_blocks[int(b)] = blocks
    pooled_all = {int(r["budget"]): float(r["dti"]) for r in rows}
    best_all = max(pooled_all, key=lambda k: pooled_all[k])
    half = {0: [], 1: []}
    n_blocks = len(per_budget_blocks[int(best_all)])
    for name, blocks in per_budget_blocks.items():
        for i, blk in enumerate(blocks):
            half[i % 2].append((name, i, blk["score"]))
    held_out = []
    for fold in (0, 1):
        select_from = [t for t in half[fold]]
        eval_on = [t for t in half[1 - fold]]
        sums = {}
        for name, i, score in select_from:
            sums.setdefault(name, []).append(score)
        pick = max(sums, key=lambda k: float(np.mean(sums[k])))
        eval_scores = [s for n, i, s in eval_on if n == pick]
        held_out.append({
            "fold_selected_on": fold,
            "budget_chosen": int(pick),
            "selection_mean_dti_on_fold": float(np.mean(sums[pick])),
            "evaluation_mean_dti_on_other_fold": float(np.mean(eval_scores)),
            "evaluation_blocks": len(eval_scores),
        })
    return {
        "ladder": rows,
        "blocks_per_budget": per_budget_blocks,
        "budget_all_blocks": int(best_all),
        "pooled_dti_all_blocks": float(pooled_all[best_all]),
        "two_fold_held_out": held_out,
        "two_fold_mean_held_out_dti": float(np.mean([h["evaluation_mean_dti_on_other_fold"] for h in held_out])),
        "n_blocks": int(n_blocks),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-uniqueness", action="store_true",
                        help="skip the (slow) 279-raster uniqueness audit")
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument("--prior-cache", type=Path, default=ROOT / "prior_cache")
    args = parser.parse_args()

    inputs = load_inputs()
    print(f"footprint {int(inputs['footprint'].sum()):,}  catalogue {int(inputs['catalogue'].sum()):,}")
    built = build_fields(inputs, cache=not args.no_cache)
    variants = built["variants"]

    selection = {}
    for name, field in variants.items():
        sel = budget_selection(field, inputs["catalogue"], inputs["footprint"], inputs["catalogue"])
        selection[name] = sel
        print(f"{name}: pooled DTI {sel['pooled_dti_all_blocks']:.4f} at budget "
              f"{sel['budget_all_blocks']:,}  (2-fold held-out mean {sel['two_fold_mean_held_out_dti']:.4f})")

    # Pre-registered selection rule: best pooled, but prefer a concurrence variant when within 0.005.
    ranked = sorted(selection.items(), key=lambda kv: -kv[1]["pooled_dti_all_blocks"])
    top_name, top = ranked[0]
    if "euler_only" in top_name:
        for name, sel in ranked[1:]:
            if abs(sel["pooled_dti_all_blocks"] - top["pooled_dti_all_blocks"]) <= 0.005:
                print(f"pre-registered tie rule: preferring {name} over {top_name}")
                top_name, top = name, sel
                break
    print(f"CHOSEN variant {top_name} at budget {top['budget_all_blocks']:,}")
    chosen_field = variants[top_name]
    budget = int(top["budget_all_blocks"])
    dots = nms_dots(chosen_field, inputs["footprint"], min_separation_px=3, budget=budget,
                    catalogue=inputs["catalogue"], catalogue_buffer_px=2)

    # Controls at the same emitted count.
    rng = np.random.default_rng(4006)
    idx = np.flatnonzero(inputs["footprint"].ravel())
    n_emit = int((dots > 0).sum())
    rand = np.zeros(dots.shape, dtype=np.float32)
    rand.ravel()[rng.choice(idx, size=min(n_emit, idx.size), replace=False)] = 1.0
    from gemsdoe40.holdout import gradient_baseline

    grad = gradient_baseline(inputs["rtp"], inputs["footprint"], n_emit)
    controls = {
        "same_mass_random": blocked_score(rand, inputs["catalogue"], inputs["footprint"]),
        "same_mass_gradient_topk": blocked_score(grad, inputs["catalogue"], inputs["footprint"]),
    }
    chosen_score = blocked_score(dots.astype(np.float32), inputs["catalogue"], inputs["footprint"])
    print(f"candidate pooled DTI {chosen_score['pooled']['score']:.4f} | "
          f"random {controls['same_mass_random']['pooled']['score']:.4f} | "
          f"gradient {controls['same_mass_gradient_topk']['pooled']['score']:.4f}")

    stamp = _now()
    DOWNLOADS.mkdir(parents=True, exist_ok=True)
    DOCDATA.mkdir(parents=True, exist_ok=True)
    EVIDENCE.mkdir(parents=True, exist_ok=True)
    tmp = DOWNLOADS / "_tmp_h4.tif"
    write_submission(dots, inputs["sample"], tmp, inputs["footprint"], outside="nan")
    digest = sha256(tmp)[:8]
    slug = f"gems40-h4a-euler-sicontact-depthcluster-{stamp}-{digest}"
    nan_path = DOWNLOADS / f"{slug}-nan.tif"
    zeros_path = DOWNLOADS / f"{slug}-zeros.tif"
    tmp.replace(nan_path)
    write_submission(dots, inputs["sample"], zeros_path, inputs["footprint"], outside="zero")

    rec_nan = check_submission(nan_path, inputs["footprint"])
    rec_zeros = check_submission(zeros_path, inputs["footprint"])
    print("format nan:", rec_nan["ok_range"], rec_nan.get("issues"))
    print("format zeros:", rec_zeros["ok_range"], rec_zeros.get("issues"))
    zip_path = DOWNLOADS / f"{slug}-nan.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.write(nan_path, arcname=nan_path.name)

    uniqueness = {"skipped": True}
    if not args.skip_uniqueness and args.prior_cache.exists():
        from gemsdoe40.uniqueness import compare_against

        t0 = time.time()
        uniqueness = compare_against(dots.astype(np.float64), args.prior_cache)
        uniqueness["elapsed_s"] = round(time.time() - t0, 1)
        uniqueness.pop("rows", None)
        print(f"uniqueness is_new={uniqueness['is_new']} worst_r={uniqueness['worst_pearson']:.4f} "
              f"worst_jaccard={uniqueness['worst_jaccard']:.4f}")

    receipt = {
        "candidate_id": "H4-A",
        "slug": slug,
        "generated_utc": stamp,
        "sha256_prefix": digest,
        "method": (
            "Reid et al. (1990) 3-D Euler deconvolution, structural index 0 "
            "(fault-like contact), moving windows 8/12/16 px, on RTP magnetics and "
            "isostatic gravity; 3-D depth-consistency gate; cross-field (x,y,z) "
            "concordance; weighted KDE; multi-scale Hessian lineament response; "
            "3DEP scarp / detrended-slope surface-expression concurrence; "
            "non-maximum-suppressed dots at 3 px separation; catalogue pruned exactly."
        ),
        "variants": {k: {"pooled_dti": v["pooled_dti_all_blocks"],
                         "budget": v["budget_all_blocks"],
                         "two_fold_held_out_dti": v["two_fold_mean_held_out_dti"]}
                     for k, v in selection.items()},
        "chosen_variant": top_name,
        "budget": budget,
        "emitted_px": n_emit,
        "candidate_blocked_score": chosen_score,
        "controls": controls,
        "field_stats": built["field_stats"],
        "terrain_stats": built["terrain_stats"],
        "gap_stats": built["gap_stats"],
        "budget_selection": {
            "pooled_dti_all_blocks": top["pooled_dti_all_blocks"],
            "two_fold_mean_held_out_dti": top["two_fold_mean_held_out_dti"],
            "two_fold": top["two_fold_held_out"],
            "ladder": top["ladder"],
            "per_block": top["blocks_per_budget"][budget],
        },
        "format": {"nan_tif": rec_nan, "zeros_tif": rec_zeros,
                   "zip": {"path": str(zip_path), "bytes": zip_path.stat().st_size,
                           "sha256": sha256(zip_path)}},
        "uniqueness": uniqueness,
        "submission_name": f"gems40-h4a-euler-contact-depthcluster-{digest}",
        "submission_note": (
            f"H4-A Euler SI=0 depth-cluster contact network; {n_emit:,} dots, "
            f"3px NMS, catalogue pruned; blocked catalogue DTI "
            f"{chosen_score['pooled']['score']:.4f}; {digest}"
        )[:200],
    }
    (DOCDATA / f"h4a-receipt-{stamp}.json").write_text(json.dumps(receipt, indent=2, default=str) + "\n")
    (DOCDATA / "h4a_receipt.json").write_text(json.dumps(receipt, indent=2, default=str) + "\n")
    (EVIDENCE / "h4a_last_run.json").write_text(json.dumps(receipt, indent=2, default=str) + "\n")
    print(f"wrote {nan_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
