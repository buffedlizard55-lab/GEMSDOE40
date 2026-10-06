#!/usr/bin/env python3
"""ARCHIVED ONLY — opt-in historical reproduction, not current submission advice.

Run dual-field Euler deconvolution → weighted KDE → legal GeoTIFF.

This is the GEMSDOE40 submission generator.  It does not copy any previous
dotted / gradient / curvature candidate.  It reads the official 19-band
feature stack, deconvolves RTP magnetics and isostatic gravity at SI=0
(fault-like contact; Reid et al. 1990), converts the depth-labelled cloud
into a weighted KDE, masks the public catalogue, writes zeros-outside and
NaN-outside GeoTIFFs, scores a 4-fold spatial holdout, and refuses to call
the file new if it correlates with a prior submission.
"""
from __future__ import annotations

import json
import sys
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe40 import BANDS, GRID_HEIGHT, GRID_WIDTH  # noqa: E402
from gemsdoe40.euler import deconvolve, merge_clouds  # noqa: E402
from gemsdoe40.grid import (  # noqa: E402
    check_submission,
    footprint_from_sample,
    read_band,
    read_labels,
    sha256,
    write_submission,
)
from gemsdoe40.holdout import evaluate, gradient_baseline  # noqa: E402
from gemsdoe40.metric import dti_binary  # noqa: E402
from gemsdoe40.kde import assemble_field  # noqa: E402
from gemsdoe40.uniqueness import compare_against  # noqa: E402

DATA = ROOT / "data"
REF = ROOT / "ref"
DOCS = ROOT / "docs"
DOWNLOADS = DOCS / "downloads"
EVIDENCE = ROOT / "evidence"


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def run() -> dict:
    from gemsdoe40.legacy_guard import require_legacy_opt_in
    require_legacy_opt_in()
    t0 = time.time()
    sample = DATA / "sample_submission.tif"
    labels = DATA / "labels.tif"
    features = DATA / "training_features.tif"
    for p in (sample, labels, features):
        if not p.exists():
            raise SystemExit(f"missing {p}")

    footprint = footprint_from_sample(sample)
    catalogue = read_labels(labels)
    assert footprint.shape == (GRID_HEIGHT, GRID_WIDTH)
    print(f"footprint {int(footprint.sum()):,}  catalogue {int(catalogue.sum()):,}")

    rtp, rtp_ok = read_band(features, BANDS["rtp"])
    grav, grav_ok = read_band(features, BANDS["iso_grav_anom"])
    # Intersection with the submission footprint (sample is the legal domain).
    rtp_ok = rtp_ok & footprint
    grav_ok = grav_ok & footprint
    print(f"rtp valid {int(rtp_ok.sum()):,}  grav valid {int(grav_ok.sum()):,}")

    clouds = []
    configs = [
        ("rtp", rtp, rtp_ok, 0.0, 8, 4),
        ("rtp", rtp, rtp_ok, 0.0, 12, 4),
        ("iso_grav_anom", grav, grav_ok, 0.0, 8, 4),
        ("iso_grav_anom", grav, grav_ok, 0.0, 12, 4),
    ]
    for name, field, valid, si, win, stride in configs:
        t1 = time.time()
        cloud = deconvolve(
            field, valid,
            field_name=name,
            structural_index=si,
            window_px=win,
            stride_px=stride,
            analytic_percentile=72.0,
            max_rel_se=0.22,
            min_depth_m=80.0,
            max_depth_m=2200.0,
        )
        print(f"  Euler {name} SI={si} win={win}: {len(cloud):,} solutions "
              f"in {time.time()-t1:.1f}s  stats={cloud.stats}")
        clouds.append(cloud)

    mag_cloud = merge_clouds([c for c in clouds if c.field_name == "rtp"])
    grav_cloud = merge_clouds([c for c in clouds if c.field_name == "iso_grav_anom"])
    print(f"merged mag {len(mag_cloud):,}  grav {len(grav_cloud):,}")

    # Skill check is on the UNMASKED KDE (catalogue still visible).  The
    # shipped file then zeros the catalogue + 2 px, because the hidden test
    # set is off-catalogue by construction.  Scoring the masked file against
    # the catalogue is a tautology (~0) and is not used as a promotion gate.
    unmasked, unmasked_stats = assemble_field(
        mag_cloud, grav_cloud, footprint, catalogue,
        catalogue_buffer_px=-1,  # no catalogue mask: skill check only
        sigma_px=1.7,
        concord_weight=1.75,
        floor_quantile=0.98,
        power=1.4,
    )
    hold = evaluate(unmasked, catalogue, footprint, grad_field=rtp, seed=40)
    print("unmasked KDE stats", json.dumps(unmasked_stats, indent=2))
    print(f"UNMASKED holdout mean blocked DTI {hold['blocked']['mean_dti']:.4f}  "
          f"full {hold['full_dti']['dti']:.4f}  "
          f"vs random Δ {hold['delta_vs_random']:.4f}  "
          f"vs gradient Δ {hold.get('delta_vs_gradient', float('nan')):.4f}")

    field, kde_stats = assemble_field(
        mag_cloud, grav_cloud, footprint, catalogue,
        catalogue_buffer_px=2,
        sigma_px=1.7,
        concord_weight=1.75,
        floor_quantile=0.98,
        power=1.4,
    )
    print("SHIPPED (catalogue-masked) KDE stats", json.dumps(kde_stats, indent=2))
    kde_stats["unmasked"] = unmasked_stats

    # Fair same-count binary contrast: top-N unmasked Euler cores vs |∇RTP|.
    n_emit = int((unmasked > 0).sum())
    u_bin = unmasked > 0
    euler_bin = dti_binary(u_bin, catalogue, valid=footprint)
    grad_bin = dti_binary(
        gradient_baseline(rtp, footprint, n_emit) > 0, catalogue, valid=footprint
    )
    hold["euler_binary_same_count"] = euler_bin
    hold["gradient_binary_same_count"] = grad_bin
    hold["delta_binary_vs_gradient"] = float(euler_bin["dti"] - grad_bin["dti"])
    print(f"binary same-count N={n_emit}  Euler DTI {euler_bin['dti']:.4f}  "
          f"grad DTI {grad_bin['dti']:.4f}  Δ {hold['delta_binary_vs_gradient']:+.4f}")

    # Uniqueness against every prior TIF sitting in ref/
    uniq = compare_against(field, REF) if REF.exists() else {"is_new": False, "n_priors": 0, "rows": [], "error": "missing prior directory; novelty unproved"}
    print(f"uniqueness is_new={uniq['is_new']}  worst_pearson={uniq.get('worst_pearson')}  "
          f"worst_jaccard={uniq.get('worst_jaccard')}  vs {uniq.get('worst_file')}")
    if not uniq["is_new"]:
        raise SystemExit(
            "REFUSING to call this new: near-duplicate of "
            f"{uniq.get('worst_file')} (pearson={uniq.get('worst_pearson')}, "
            f"jaccard={uniq.get('worst_jaccard')})"
        )

    stamp = _now()
    DOWNLOADS.mkdir(parents=True, exist_ok=True)
    EVIDENCE.mkdir(parents=True, exist_ok=True)

    # Write zeros (portal-safe) and nan (sample-format) twins.
    tmp_zeros = DOWNLOADS / "_tmp_zeros.tif"
    write_submission(field, sample, tmp_zeros, footprint, outside="zero")
    content8 = sha256(tmp_zeros)[:8]
    slug = f"gemsdoe40-euler-si0-depthkde-{stamp}-{content8}"
    zeros_path = DOWNLOADS / f"{slug}-zeros.tif"
    nan_path = DOWNLOADS / f"{slug}-nan.tif"
    tmp_zeros.replace(zeros_path)
    write_submission(field, sample, nan_path, footprint, outside="nan")

    rec_z = check_submission(zeros_path, footprint)
    rec_n = check_submission(nan_path, footprint)
    # Portal guarantee: every cell of the zeros file is finite and in [0, 1].
    with __import__("rasterio").open(zeros_path) as s:
        a = s.read(1)
    assert a.shape == (GRID_HEIGHT, GRID_WIDTH)
    assert a.dtype == np.float32
    assert np.isfinite(a).all()
    assert float(a.min()) >= 0.0 and float(a.max()) <= 1.0
    assert s.nodata is None
    assert rec_z["ok_range"] and rec_n["ok_range"]
    assert kde_stats["on_catalogue_positive"] == 0

    zip_path = DOWNLOADS / f"{slug}-zeros.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.write(zeros_path, arcname=zeros_path.name)

    note = (
        f"GEMSDOE40 Euler SI=0 depth-KDE | RTP+isograv win8/12, "
        f"{kde_stats['emitted_px']} px, 0 on cat, "
        f"bin Δvs|∇RTP| {hold.get('delta_binary_vs_gradient', 0):+.4f} | {content8}"
    )
    if len(note) > 200:
        note = note[:197] + "..."

    receipt = {
        "candidate_id": "H40-EULER-SI0-DEPTHKDE",
        "slug": slug,
        "content_digest8": content8,
        "method": (
            "Reid et al. 1990 3-D Euler deconvolution, SI=0 (fault-like contact), "
            "on RTP magnetics and isostatic gravity; weighted KDE of shallow, "
            "tight, mutually-consistent solution clusters; catalogue-buffered to 0."
        ),
        "citation": {
            "reid1990": "https://doi.org/10.1190/1.1442774",
            "reid1990_pdf": "https://www.reid-geophys.co.uk/wp-content/uploads/2017/11/Reid-et-al-1990.pdf",
            "problem": "https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/",
        },
        "submission_name": "GEMSDOE40-euler-si0-depthkde",
        "submission_note": note,
        "zeros_tif": rec_z,
        "nan_tif": rec_n,
        "zip": {"path": str(zip_path), "bytes": zip_path.stat().st_size, "sha256": sha256(zip_path)},
        "kde": kde_stats,
        "euler_mag": mag_cloud.stats,
        "euler_grav": grav_cloud.stats,
        "holdout": {
            "mean_blocked_dti": hold["blocked"]["mean_dti"],
            "full_dti": hold["full_dti"],
            "random_dti": hold["random_same_mass_dti"],
            "delta_vs_random": hold["delta_vs_random"],
            "gradient_dti": hold.get("gradient_same_mass_dti"),
            "delta_vs_gradient": hold.get("delta_vs_gradient"),
            "per_fold": {k: v for k, v in hold["blocked"]["per_fold"].items()},
        },
        "uniqueness": uniq,
        "elapsed_s": round(time.time() - t0, 1),
        "generated_utc": stamp,
    }
    def _jsonable(obj):
        if isinstance(obj, dict):
            return {k: _jsonable(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [_jsonable(v) for v in obj]
        if isinstance(obj, float) and (obj != obj):  # NaN
            return None
        if hasattr(obj, "item"):
            try:
                return _jsonable(obj.item())
            except Exception:
                return str(obj)
        return obj

    rec_path = DOWNLOADS / f"{slug}-audit.json"
    payload = _jsonable(receipt)
    rec_path.write_text(json.dumps(payload, indent=2))
    (EVIDENCE / "last_run.json").write_text(json.dumps(payload, indent=2))
    print(f"wrote {zeros_path}  sha256 {rec_z['sha256']}")
    print(f"note ({len(note)}/200): {note}")
    print(f"elapsed {receipt['elapsed_s']}s")
    return receipt


if __name__ == "__main__":
    run()
