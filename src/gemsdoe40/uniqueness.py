"""Hash and correlate a candidate raster against every prior submission.

The brief requires that we refuse to call a file new if it is a near-duplicate
of a gradient / curvature / dotted-H19 candidate already published.  Depth
clustering should produce a *visibly different spatial pattern*.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio

NEAR_DUPLICATE_CORR = 0.85
NEAR_DUPLICATE_JACCARD = 0.50


def _load(path: str | Path) -> np.ndarray:
    with rasterio.open(path) as s:
        a = s.read(1).astype(np.float64)
    a = np.where(np.isfinite(a), a, 0.0)
    return a


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    x = a.ravel()
    y = b.ravel()
    if x.size != y.size:
        raise ValueError("shape mismatch")
    x = x - x.mean()
    y = y - y.mean()
    d = float(np.linalg.norm(x) * np.linalg.norm(y))
    if d == 0:
        return 0.0
    return float(np.dot(x, y) / d)


def jaccard_positive(a: np.ndarray, b: np.ndarray) -> float:
    pa = a > 0
    pb = b > 0
    inter = int((pa & pb).sum())
    union = int((pa | pb).sum())
    return float(inter / union) if union else 0.0


def compare_against(candidate: np.ndarray, prior_dir: str | Path) -> dict:
    """Correlate ``candidate`` with every GeoTIFF in ``prior_dir``.

    Returns a report.  ``is_new`` is False if any prior file exceeds the
    near-duplicate thresholds.
    """
    prior_dir = Path(prior_dir)
    files = sorted(prior_dir.glob("*.tif"))
    rows = []
    worst_corr = 0.0
    worst_jac = 0.0
    worst_name = None
    for f in files:
        try:
            b = _load(f)
        except Exception as exc:  # pragma: no cover
            rows.append({"file": f.name, "error": str(exc)})
            continue
        if b.shape != candidate.shape:
            rows.append({"file": f.name, "error": f"shape {b.shape}"})
            continue
        c = pearson(candidate, b)
        j = jaccard_positive(candidate, b)
        rows.append({
            "file": f.name,
            "pearson": round(c, 6),
            "jaccard_positive": round(j, 6),
            "prior_positive_px": int((b > 0).sum()),
        })
        if abs(c) >= abs(worst_corr):
            worst_corr, worst_name = c, f.name
        if j >= worst_jac:
            worst_jac = j
    is_new = (abs(worst_corr) < NEAR_DUPLICATE_CORR) and (worst_jac < NEAR_DUPLICATE_JACCARD)
    return {
        "n_priors": len(files),
        "is_new": bool(is_new),
        "worst_pearson": float(worst_corr),
        "worst_jaccard": float(worst_jac),
        "worst_file": worst_name,
        "threshold_pearson": NEAR_DUPLICATE_CORR,
        "threshold_jaccard": NEAR_DUPLICATE_JACCARD,
        "rows": rows,
    }
