#!/usr/bin/env python3
"""H11 pre-screen: seismicity-corridor x shallow-Euler intersection on the frozen proxy.

Implements the frozen protocol in docs/research/h8-preregistration-20261006.md,
"2026-10-06 addendum 4 — H11 pre-screen". Measurement only: never writes a submission
candidate, never changes a gate. Fails closed if the preregistration document drifted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import rasterio
from scipy.ndimage import gaussian_filter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from gemsdoe40.raster import band_index_by_name  # noqa: E402
from gemsdoe40.research_holdout import read_proxy_truth, score_array_on_proxy  # noqa: E402
from gemsdoe40.trace_lock import CLOUD_SHA256, load_cloud_csv  # noqa: E402
from audit_trace_candidate import mass_matched_controls  # noqa: E402

PREREG_PATH = ROOT / "docs/research/h8-preregistration-20261006.md"
PINNED_PREREG_SHA256 = "d6117108c59f18371f25a428b608193d98abc0b03fd13f1c018d16b9c76a2a99"
FEATURES_SHA256 = "4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5"
TEMPLATE_SHA256 = "2176d08e485aa2cd2860ce8df539db4faf4d76163b38a4dd8c30a40454d35cbc"
CLOUD_PATH = ROOT / "docs/downloads/h4-euler-solutions.csv.gz"
BANDS = ("deq_n100a15", "ieq_n100a15")
SIGMA_PX = 2.0  # frozen: 200 m Gaussian before derivatives / KDE kernel
BAR_ABSOLUTE = 0.039354  # frozen: 2 x H13 session-2 G3 score 0.019677
BAR_VS_RANDOM = 2.0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def percentile_normalize(field: np.ndarray) -> np.ndarray:
    out = np.full(field.shape, np.nan, dtype=np.float64)
    finite = np.isfinite(field)
    vals = field[finite]
    if vals.size < 2:
        out[finite] = 0.0
        return out
    order = np.argsort(vals, kind="stable")
    ranks = np.empty_like(order)
    ranks[order] = np.arange(vals.size)
    out[finite] = ranks / (vals.size - 1)
    return out


def horizontal_gradient(field: np.ndarray) -> np.ndarray:
    gy, gx = np.gradient(field)
    return np.hypot(gx, gy)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, default=ROOT / "work" / "h11")
    args = parser.parse_args()

    prereg_sha = sha256(PREREG_PATH)
    if prereg_sha != PINNED_PREREG_SHA256:
        print(f"ABORT: preregistration document changed ({prereg_sha[:12]}... != pinned)")
        return 3

    features = ROOT / "data" / "training_features.tif"
    if sha256(features) != FEATURES_SHA256:
        print("ABORT: training_features.tif sha mismatch")
        return 3
    template = ROOT / "data" / "sample_submission.tif"
    if sha256(template) != TEMPLATE_SHA256:
        print("ABORT: sample_submission.tif sha mismatch")
        return 3
    if sha256(CLOUD_PATH) != CLOUD_SHA256:
        print("ABORT: frozen Euler cloud sha mismatch")
        return 3

    started = datetime.now(timezone.utc).isoformat()
    with rasterio.open(template) as tds:
        template_band = tds.read(1)
        transform = tds.transform
    footprint_mask = np.isfinite(template_band)

    with rasterio.open(features) as ds:
        raw: dict[str, np.ndarray] = {}
        for name in BANDS:
            band = ds.read(band_index_by_name(ds, name)).astype(np.float64)
            band[~np.isfinite(band)] = np.nan
            raw[name] = band

    smooth = {n: gaussian_filter(np.nan_to_num(b, nan=0.0), SIGMA_PX, mode="nearest")
              for n, b in raw.items()}
    for n in smooth:
        smooth[n][~footprint_mask] = np.nan

    # Frozen Euler-support KDE from the hash-pinned H4 cloud (weighted records only),
    # using the exact H8-family construction: bilinear splat + Gaussian sigma=2 truncate=4.
    from gemsdoe40.contact_euler import bilinear_splat
    cloud = load_cloud_csv(CLOUD_PATH, transform)  # structured: family,row,col,depth_m,weight,...
    rows = cloud["row"]
    cols = cloud["col"]
    weights = cloud["weight"]
    support = bilinear_splat(rows, cols, weights, footprint_mask.shape)
    euler_kde = gaussian_filter(support, sigma=2.0, truncate=4.0, mode="constant", cval=0.0)
    euler_kde[~footprint_mask] = np.nan

    prox_deq = percentile_normalize(-raw["deq_n100a15"])
    prox_ieq = percentile_normalize(-raw["ieq_n100a15"])
    hg_deq = horizontal_gradient(smooth["deq_n100a15"])
    hg_ieq = horizontal_gradient(smooth["ieq_n100a15"])
    prox_both = prox_deq * prox_ieq
    rank_euler = percentile_normalize(euler_kde)
    prox_x_euler = prox_deq * rank_euler

    indicators = {
        "prox_deq": prox_deq,
        "prox_ieq": prox_ieq,
        "HG_deq": hg_deq,
        "HG_ieq": hg_ieq,
        "prox_both": prox_both,
        "prox_deq_x_euler": prox_x_euler,
    }

    truth, valid, labels, proxy_info = read_proxy_truth(
        ROOT / "data/external/derived_sgmc_faults_100m_u8.tif",
        template,
        ROOT / "data/labels.tif",
    )

    results: dict[str, dict] = {}
    with rasterio.open(features) as ds:
        tmi = ds.read(band_index_by_name(ds, "tmi")).astype(np.float64)
    tmi[~np.isfinite(tmi)] = 0.0

    for name, field in indicators.items():
        emission = percentile_normalize(field).astype(np.float32)
        score = score_array_on_proxy(emission, truth, valid, labels)["pooled"]["score"]
        controls = mass_matched_controls(emission, footprint_mask, tmi)
        control_scores = {cname: score_array_on_proxy(cfield, truth, valid, labels)["pooled"]["score"]
                          for cname, cfield in controls["fields"].items()}
        best_random = max(v for k, v in control_scores.items() if k.startswith("random_mass_matched"))
        results[name] = {
            "positive_px": int((emission > 0).sum()),
            "mass_sum": float(np.nansum(emission)),
            "mass_budget_px": controls["mass_budget_px"],
            "proxy_score": float(score),
            "control_scores": {k: float(v) for k, v in control_scores.items()},
            "best_random_control": float(best_random),
            "passes_absolute_bar": bool(score > BAR_ABSOLUTE),
            "passes_random_bar": bool(score > BAR_VS_RANDOM * best_random),
        }

    best = max(results, key=lambda k: results[k]["proxy_score"])
    advance = results[best]["passes_absolute_bar"] and results[best]["passes_random_bar"]

    out = {
        "protocol": "h11-prescore frozen addendum 4 (2026-10-06)",
        "preregistration_sha256_at_runtime": prereg_sha,
        "script_sha256": sha256(Path(__file__).resolve()),
        "features_sha256": sha256(features),
        "cloud_sha256": sha256(CLOUD_PATH),
        "proxy_info": proxy_info,
        "started_utc": started,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "sigma_px": SIGMA_PX,
        "bars": {"absolute": BAR_ABSOLUTE, "vs_random_factor": BAR_VS_RANDOM},
        "indicators": results,
        "best_indicator": best,
        "decision": "ADVANCE-TO-FULL-H11" if advance else "NEGATIVE-AT-PRESCREEN",
        "slot_used": False,
    }
    args.work.mkdir(parents=True, exist_ok=True)
    target = args.work / "prescore.json"
    target.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"decision": out["decision"], "best": best,
                      "best_score": results[best]["proxy_score"],
                      "best_random": results[best]["best_random_control"]}, indent=2))
    print(f"wrote {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
