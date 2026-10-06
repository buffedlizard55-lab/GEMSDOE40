#!/usr/bin/env python3
"""Re-measure the live-score calibration instrument on the anchor priors' own bytes.

The sibling repository's saturating model

    TP_hat = K (1 - exp(-a n w / K))
    S_hat  = TP_hat / (0.2 TP_hat + 0.2 n + 0.8 K)

was fitted to 12 organizer-scored artifacts with ``w`` = mean per-dot kernel credit
against *its own* surrogate definition.  This repository's instrument audit found the
form does not rank the scored artifacts when ``w`` is re-measured on the off-catalogue
SGMC surrogate.  This script does the measurement itself:

* every anchor raster in ``ref/prior`` is byte-checked against the pinned inventory,
* ``n`` (emitted mass off-catalogue, inside the scoreable footprint) and ``w`` are
  measured from the pixels,
* the same functional form is **refitted** by non-linear least squares on the measured
  pairs and scored by leave-one-out Spearman and RMSE,
* the inherited constants are evaluated on the same data for comparison.

Nothing here is a competition score.  Output: ``docs/data/instrument-recalibration.json``.
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
from scipy import ndimage, optimize, stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe40.emission import INSTRUMENT_A, INSTRUMENT_K, calibrated_live_score  # noqa: E402
from gemsdoe40.measure import credit_components  # noqa: E402

DATA = ROOT / "data"
PRIOR = ROOT / "ref" / "prior"
CORPUS = ROOT / "ref" / "meta" / "prior_corpus.md"
OUT = ROOT / "docs" / "data" / "instrument-recalibration.json"
RADIUS_PX = 3.0


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_corpus() -> dict[str, dict]:
    """short name -> {sha16, mass, w, live, pred} from the committed measured table."""
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
            "corpus_mass": int(m.group(4).replace(",", "")),
            "corpus_w": float(m.group(5)),
            "live_score": None if live_raw in ("—", "-", "") else float(live_raw),
            "corpus_pred": float(m.group(7)),
        }
    return rows


def load_masks() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with rasterio.open(DATA / "sample_submission.tif") as ds:
        footprint = np.isfinite(ds.read(1))
    with rasterio.open(DATA / "labels.tif") as ds:
        catalogue = ds.read(1) > 0
    with rasterio.open(DATA / "external" / "derived_sgmc_faults_100m_u8.tif") as ds:
        sgmc = ds.read(1) > 0
    return footprint, catalogue, sgmc


def measure_dots(scored: np.ndarray, truth: np.ndarray, truth_kernel: np.ndarray) -> dict:
    """Measure the instrument's inputs from the prediction's own pixels."""
    dt = ndimage.distance_transform_edt(~truth)
    k_truth = np.maximum(1.0 - dt / RADIUS_PX, 0.0)
    del dt
    w = float(k_truth[scored].mean())
    credit_pred = float((truth_kernel * scored).sum())
    # ``u`` is the mirror quantity: credit the prediction delivers to the surrogate truth.
    return {"w": w, "u": credit_pred / float(scored.sum()), "credit_pred": credit_pred}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()

    corpus = parse_corpus()
    footprint, catalogue, sgmc = load_masks()
    truth = sgmc & ~catalogue & footprint
    print(f"footprint {int(footprint.sum()):,}  catalogue {int(catalogue.sum()):,}  "
          f"SGMC {int(sgmc.sum()):,}  surrogate truth (SGMC off-catalogue) {int(truth.sum()):,}")

    dt_truth = ndimage.distance_transform_edt(~truth)
    truth_kernel = np.maximum(1.0 - dt_truth / RADIUS_PX, 0.0).astype(np.float32)
    del dt_truth

    by_sha: dict[str, Path] = {}
    for p in sorted(PRIOR.glob("*.tif")):
        by_sha[sha256(p)[:16]] = p

    rows: list[dict] = []
    for short, meta in corpus.items():
        path = by_sha.get(meta["sha16"])
        if path is None:
            continue
        with rasterio.open(path) as ds:
            a = ds.read(1)
        inside = np.isfinite(a) & (a > 0) & footprint
        scored = inside & ~catalogue
        n = int(scored.sum())
        if n:
            m = measure_dots(scored, truth, truth_kernel)
        else:
            m = {"w": 0.0, "u": 0.0, "credit_pred": 0.0}
        comp = credit_components(scored.astype(np.float32), truth, valid=footprint)
        rows.append({
            "proxy_dti": comp["dti"], "proxy_tp": comp["tp"], "proxy_fp": comp["fp"],
            "proxy_fn": comp["fn"], "proxy_truth_px": comp["n_truth"],
            "credit_per_mass_u": comp["credit_per_mass"],
            "anchor": short, "file": path.name, "sha256": sha256(path),
            "bytes": path.stat().st_size,
            "mass": n, "on_catalogue_px": int((inside & catalogue).sum()),
            "corpus_mass": meta["corpus_mass"], "corpus_w": meta["corpus_w"],
            "w": m["w"], "u": m["u"], "credit_pred": m["credit_pred"],
            "live_score": meta["live_score"], "corpus_pred": meta["corpus_pred"],
        })
        ratio = (meta["live_score"] / comp["dti"]) if meta["live_score"] and comp["dti"] > 0 else None
        print(f"{short:52s} n={n:7,d} w={m['w']:.4f} proxyDTI={comp['dti']:.4f} "
              f"live={meta['live_score']} ratio={None if ratio is None else round(ratio, 2)}")
        del a, inside, scored

    scored_rows = [r for r in rows if r["live_score"] is not None]
    n_all, w_all = np.array([r["mass"] for r in scored_rows], float), \
        np.array([r["w"] for r in scored_rows], float)
    y = np.array([r["live_score"] for r in scored_rows], float)

    def model(params, n, w):
        K, a = np.exp(params)  # positive by construction
        tp = K * (1.0 - np.exp(-a * n * w / K))
        return tp / (0.2 * tp + 0.2 * n + 0.8 * K)

    def sse(params):
        return float(np.sum((model(params, n_all, w_all) - y) ** 2))

    fitted = optimize.minimize(sse, np.log([INSTRUMENT_K, INSTRUMENT_A]), method="Nelder-Mead",
                               options=dict(xatol=1e-8, fatol=1e-12, maxiter=20000))
    K_hat, a_hat = np.exp(fitted.x)
    pred_in = model(fitted.x, n_all, w_all)

    # leave-one-out refits
    loo_pred, loo_params = [], []
    for i in range(len(scored_rows)):
        keep = np.arange(len(scored_rows)) != i
        res = optimize.minimize(
            lambda p: float(np.sum((model(p, n_all[keep], w_all[keep]) - y[keep]) ** 2)),
            np.log([INSTRUMENT_K, INSTRUMENT_A]), method="Nelder-Mead",
            options=dict(xatol=1e-8, fatol=1e-12, maxiter=20000))
        loo_params.append(np.exp(res.x).tolist())
        loo_pred.append(float(model(res.x, n_all[i:i + 1], w_all[i:i + 1])[0]))

    inherited = np.array([calibrated_live_score(r["mass"], r["w"]) for r in scored_rows])
    loo_pred_arr = np.array(loo_pred)

    def spear(a_, b_):
        if len(a_) < 3 or np.ptp(a_) == 0 or np.ptp(b_) == 0:
            return None
        return float(stats.spearmanr(a_, b_).statistic)

    # ---- live/proxy transfer calibration ---------------------------------- #
    proxy = np.array([r["proxy_dti"] for r in scored_rows], float)
    transfer = {}
    for name, transform in (("linear_origin", "lin"), ("loglog", "log")):
        if name == "loglog":
            if (proxy <= 0).any() or (y <= 0).any():
                continue
            slope, intercept = np.polyfit(np.log(proxy), np.log(y), 1)
            fit_fn = lambda p_: float(np.exp(intercept) * p_ ** slope)  # noqa: E731
            params = {"a": float(np.exp(intercept)), "b": float(slope)}
        else:
            a = float((proxy @ y) / (proxy @ proxy))
            fit_fn = lambda p_: float(a * p_)  # noqa: E731
            params = {"a": a}
        preds_in = np.array([fit_fn(p_) for p_ in proxy])
        loo = []
        for i in range(len(proxy)):
            keep = np.arange(len(proxy)) != i
            if name == "loglog":
                sl, ic = np.polyfit(np.log(proxy[keep]), np.log(y[keep]), 1)
                loo.append(float(np.exp(ic) * proxy[i] ** sl))
            else:
                aa = float((proxy[keep] @ y[keep]) / (proxy[keep] @ proxy[keep]))
                loo.append(float(aa * proxy[i]))
        transfer[name] = {
            **params,
            "in_sample_spearman": spear(preds_in, y),
            "in_sample_rmse": float(np.sqrt(np.mean((preds_in - y) ** 2))),
            "max_abs_error": float(np.max(np.abs(preds_in - y))),
            "loo_spearman": spear(np.array(loo), y),
            "loo_rmse": float(np.sqrt(np.mean((np.array(loo) - y) ** 2))),
        }
    # ---- hidden-truth transfer model -------------------------------------- #
    # Exact metric algebra with two unknowns: the hidden truth has N pixels, and the
    # credit a dot earns against it is lambda times the credit it earns against the
    # surrogate (u_hidden = lambda * w_surrogate, and C = M*w by definition of w).
    #     DTI = lambda*M*w / (0.2*lambda*M*w + 0.2*M*(1 - lambda*w) + 0.8*N)
    def transfer_model(params, mass, w):
        lam, n_true = np.exp(params)
        t = lam * mass * w
        c = np.minimum(lam * w, 1.0) * mass
        den = 0.2 * t + 0.2 * (mass - c) + 0.8 * n_true
        return np.where(den > 0, t / np.maximum(den, 1e-12), 0.0)

    from scipy import optimize as _opt

    def _sse(p):
        return float(np.sum((transfer_model(p, n_all, w_all) - y) ** 2))

    hidden_fit = _opt.minimize(_sse, np.log([1.0, 1e4]), method="Nelder-Mead",
                               options=dict(xatol=1e-10, fatol=1e-14, maxiter=40000))
    lam_hat, n_hat = np.exp(hidden_fit.x)
    pred_h = transfer_model(hidden_fit.x, n_all, w_all)
    loo_h = []
    for i in range(len(scored_rows)):
        keep = np.arange(len(scored_rows)) != i
        res = _opt.minimize(
            lambda pp: float(np.sum((transfer_model(pp, n_all[keep], w_all[keep]) - y[keep]) ** 2)),
            np.log([1.0, 1e4]), method="Nelder-Mead",
            options=dict(xatol=1e-10, fatol=1e-14, maxiter=40000))
        loo_h.append(float(transfer_model(res.x, n_all[i:i + 1], w_all[i:i + 1])[0]))
    hidden_transfer = {
        "lambda": float(lam_hat), "N_hidden_px": float(n_hat),
        "sse": float(hidden_fit.fun),
        "in_sample_spearman": spear(pred_h, y),
        "in_sample_rmse": float(np.sqrt(np.mean((pred_h - y) ** 2))),
        "max_abs_error": float(np.max(np.abs(pred_h - y))),
        "loo_spearman": spear(np.array(loo_h), y),
        "loo_rmse": float(np.sqrt(np.mean((np.array(loo_h) - y) ** 2))),
        "loo_pred": loo_h,
        "form": ("DTI = lam*M*w / (0.2*lam*M*w + 0.2*M*(1 - min(lam*w,1)) + 0.8*N); "
                 "exact metric algebra with two fitted unknowns"),
    }
    # Per-anchor instrument readings under the refit, the inherited constants and the
    # two transfer models — the table the site quotes must come from here, not from prose.
    instrument_table = {}
    for r in rows:  # scored anchors plus the unscored comparators, same code path
        instrument_table[r["anchor"]] = {
            "mass": r["mass"], "w": r["w"], "proxy_dti": r["proxy_dti"],
            "live_score": r.get("live_score"),
            "instrument": float(calibrated_live_score(r["mass"], r["w"], K=K_hat, a=a_hat)),
            "instrument_inherited": float(calibrated_live_score(r["mass"], r["w"])),
            "ratio_live_over_proxy": (r["live_score"] / r["proxy_dti"]
                                      if r.get("live_score") and r["proxy_dti"] > 0 else None),
        }
    ratios = [r["live_score"] / r["proxy_dti"] for r in scored_rows if r["proxy_dti"] > 0]
    report = {
        "live_over_proxy_dti_ratio": {
            "mean": float(np.mean(ratios)), "median": float(np.median(ratios)),
            "min": float(np.min(ratios)), "max": float(np.max(ratios)),
            "note": ("empirical transfer of the off-catalogue SGMC surrogate DTI to the "
                     "owner-reported live scores; not a causal relation and not evidence "
                     "about the organizer's hidden set"),
        },
        "transfer_calibration": transfer,
        "instrument_table": instrument_table,
        "hidden_truth_transfer_model": hidden_transfer,
        "generated_from": "ref/prior anchor bytes + data/sample_submission.tif + data/labels.tif + "
                          "data/external/derived_sgmc_faults_100m_u8.tif",
        "surrogate_truth_px": int(truth.sum()),
        "n_anchors_measured": len(rows),
        "n_anchors_with_live_score": len(scored_rows),
        "score_provenance": ("live values are the owner-reported associations recorded in "
                             "ref/meta/prior_corpus.md; they are not organizer-authenticated here"),
        "measurement": {
            "mass": "emitted pixels with value > 0 inside the sample footprint and off-catalogue",
            "w": "mean over emitted dots of max(0, 1 - d/3px) to the nearest SGMC off-catalogue pixel",
            "u": "truth-side credit per emitted dot on the same surrogate",
        },
        "refit": {
            "K": float(K_hat), "a": float(a_hat),
            "sse": float(fitted.fun),
            "in_sample_spearman": spear(pred_in, y),
            "in_sample_rmse": float(np.sqrt(np.mean((pred_in - y) ** 2))),
            "max_abs_error": float(np.max(np.abs(pred_in - y))),
            "loo_spearman": spear(loo_pred_arr, y),
            "loo_rmse": float(np.sqrt(np.mean((loo_pred_arr - y) ** 2))),
            "converged": bool(fitted.success),
        },
        "inherited_constants": {
            "K": INSTRUMENT_K, "a": INSTRUMENT_A,
            "in_sample_spearman": spear(inherited, y),
            "in_sample_rmse": float(np.sqrt(np.mean((inherited - y) ** 2))),
            "max_abs_error": float(np.max(np.abs(inherited - y))),
        },
        "rows": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in report.items() if k != "rows"}, indent=2))
    print(f"[out] {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
