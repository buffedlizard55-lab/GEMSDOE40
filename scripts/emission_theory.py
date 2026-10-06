#!/usr/bin/env python3
"""Metric-optimal emission study for the DOE GEMS distance-weighted Tversky index.

Question answered here, from the official metric definition and from the
organizer-supplied catalogue raster on the exact competition grid:

    For a *given* trace hypothesis, how should prediction mass be emitted so
    that the distance-weighted Tversky index is maximised?

Official metric (verbatim structure, page 967):

    DTI = TPw / (TPw + alpha*FPw + beta*FNw),   alpha=0.2, beta=0.8, R=300 m
    TPw = sum_{g in G} max_{x <= R} p(x) k(d(x,g))
    FPw = sum_{x: p(x)>0} p(x) [1 - max_{g} k(d(x,g))]
    FNw = |G| - TPw

Two consequences that are exact, not approximate:

  1.  TPw + alpha*FPw + beta*FNw = alpha*(TPw+FPw) + beta*|G|
      (substitute FNw = |G| - TPw and use alpha+beta = 1).
      With alpha=0.2, beta=0.8 this is 0.2*(TPw+FPw) + 0.8*|G|.
      => the DTI is a function of exactly two predicted quantities:
         (TPw+FPw)  and  TPw.

  2.  A unit of prediction mass placed on a trace can be credited to more
      than one ground-truth cell: one dot on a 1-px trace is the arg-max for
      every truth cell within R, and its credit is sum_d k(d).  Mass placed
      every pixel along the trace is credited once per truth cell.  Therefore
      *sparse* mass along a known trace is far more efficient than dense mass.

This script measures consequence 2 on real fault geometry: the organizer's
own catalogue raster (60,988 positive 100 m cells inside the sample footprint)
is used as the trace geometry, dots are placed every `s` pixels along those
traces, and the exact official metric is evaluated.  A jitter experiment
offsets the whole dot set to emulate real locational error.

Outputs `docs/data/emission_theory.json`.
"""
from __future__ import annotations

import json
import sys
from collections import deque
from pathlib import Path

import numpy as np
import rasterio
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe40.metric import dti_binary  # noqa: E402

TRACE_SPACINGS = (1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 24)
JITTERS_PX = (0.0, 1.0, 2.0)


def load(path: Path) -> np.ndarray:
    with rasterio.open(path) as src:
        return src.read(1)


def path_order(mask: np.ndarray, start: tuple[int, int]) -> list[tuple[int, int]]:
    """Depth-first walk of a thin component; returns pixel order along the trace."""
    h, w = mask.shape
    seen = np.zeros_like(mask, dtype=bool)
    order: list[tuple[int, int]] = []
    stack = [start]
    while stack:
        r, c = stack.pop()
        if seen[r, c]:
            continue
        seen[r, c] = True
        order.append((r, c))
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                if dr == 0 and dc == 0:
                    continue
                nr, nc = r + dr, c + dc
                if 0 <= nr < h and 0 <= nc < w and mask[nr, nc] and not seen[nr, nc]:
                    stack.append((nr, nc))
    return order


def endpoints(mask: np.ndarray) -> list[tuple[int, int]]:
    """8-neighbour count == 1 identifies trace endpoints."""
    kernel = np.ones((3, 3), dtype=np.uint8)
    kernel[1, 1] = 0
    neigh = ndimage.convolve(mask.astype(np.uint8), kernel, mode="constant", cval=0)
    ys, xs = np.nonzero(mask & (neigh <= 1))
    return list(zip(ys.tolist(), xs.tolist()))


def dotted_traces(traces: np.ndarray, spacing: int,
                  rng: np.random.Generator | None = None, jitter_px: float = 0.0) -> np.ndarray:
    """Sample every ``spacing``-th pixel along each trace component."""
    out = np.zeros(traces.shape, dtype=bool)
    labels, n = ndimage.label(traces, structure=np.ones((3, 3)))
    starts = endpoints(traces)
    by_label: dict[int, tuple[int, int]] = {}
    for r, c in starts:
        by_label.setdefault(int(labels[r, c]), (r, c))
    for lab in range(1, n + 1):
        comp = labels == lab
        start = by_label.get(lab)
        if start is None:  # closed loop: take any pixel
            rr, cc = np.nonzero(comp)
            start = (int(rr[0]), int(cc[0]))
        order = path_order(comp, start)
        keep = order[::spacing]
        for r, c in keep:
            if jitter_px <= 0:
                out[r, c] = True
            else:
                dr = int(round(rng.normal(0.0, jitter_px)))
                dc = int(round(rng.normal(0.0, jitter_px)))
                rr, cc = r + dr, c + dc
                if 0 <= rr < out.shape[0] and 0 <= cc < out.shape[1]:
                    out[rr, cc] = True
    return out


def main() -> int:
    labels = load(ROOT / "data" / "labels.tif")
    sample = load(ROOT / "data" / "sample_submission.tif")
    footprint = np.isfinite(sample)
    traces = (labels == 1) & footprint

    labels_lbl, n_comp = ndimage.label(traces, structure=np.ones((3, 3)))
    sizes = np.bincount(labels_lbl.ravel())[1:]
    order = path_order(traces, endpoints(traces)[0]) if traces.any() else []
    del order

    baseline = dti_binary(traces, traces, valid=footprint)
    result: dict = {
        "metric": "distance-weighted Tversky index, alpha=0.2, beta=0.8, R=3 px (300 m)",
        "identity": "DTI = TPw / (0.2*(TPw+FPw) + 0.8*|G|)",
        "truth_used": {
            "source": "data/labels.tif (organizer existing-fault catalogue, float32 twin = data/sample_submission.tif)",
            "positive_cells": int(traces.sum()),
            "trace_components_8conn": int(n_comp),
            "median_component_px": float(np.median(sizes)),
            "footprint_cells": int(footprint.sum()),
        },
        "perfect_dense_trace_prediction": baseline,
        "spacing_sweep": [],
    }

    rng = np.random.default_rng(20261006)
    for jitter in JITTERS_PX:
        for spacing in TRACE_SPACINGS:
            dots = dotted_traces(traces, spacing,
                                 rng=rng if jitter > 0 else None, jitter_px=jitter)
            res = dti_binary(dots, traces, valid=footprint)
            result["spacing_sweep"].append({
                "jitter_sigma_px": jitter,
                "dot_spacing_px": spacing,
                "emitted_px": int(dots.sum()),
                "tpw": res["tp"],
                "fpw": res["fp"],
                "fnw": res["fn"],
                "dti": res["dti"],
                "credit_per_dot": float(res["tp"] / max(int(dots.sum()), 1)),
                "credit_per_unit_mass": float(res["tp"] / (0.2 * (res["tp"] + res["fp"]))),
            })
            print(f"jitter={jitter:>3}  spacing={spacing:>2}  dots={int(dots.sum()):>7,}  "
                  f"TPw={res['tp']:>9.0f}  FPw={res['fp']:>8.0f}  DTI={res['dti']:.4f}  "
                  f"credit/dot={res['tp']/max(int(dots.sum()),1):.2f}")

    out = ROOT / "docs" / "data" / "emission_theory.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {out}")

    # Nearest-neighbour spacing of the owner-reported 0.2778 file, if present.
    prior = ROOT / "work" / "priors" / "h33-2-b2-nan.tif"
    if prior.exists():
        with rasterio.open(prior) as src:
            arr = src.read(1)
        pos = np.isfinite(arr) & (arr > 0)
        ys, xs = np.nonzero(pos)
        from scipy.spatial import cKDTree
        tree = cKDTree(np.column_stack([ys, xs]))
        d, _ = tree.query(np.column_stack([ys, xs]), k=2)
        nn = d[:, 1]
        result["prior_h33_2_b2"] = {
            "path": str(prior.relative_to(ROOT)),
            "positive_cells": int(pos.sum()),
            "nn_distance_px_mean": float(nn.mean()),
            "nn_distance_px_median": float(np.median(nn)),
            "nn_distance_px_p05": float(np.percentile(nn, 5)),
            "nn_distance_px_p95": float(np.percentile(nn, 95)),
            "fraction_nn_below_2px": float((nn < 2).mean()),
            "isolated_components": int(ndimage.label(pos, structure=np.ones((3, 3)))[1]),
        }
        print("h33-2-b2 nn spacing:", json.dumps(result["prior_h33_2_b2"], indent=1))
        out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
