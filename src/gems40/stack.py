"""Construction of the cross-family Euler depth-cluster field.

Kept in ``src/`` (rather than in the build script) because two callers must agree **exactly**:
``scripts/build_submission.py`` (what gets shipped) and ``scripts/multiscale_test.py`` (what gets
measured).  Duplicating the formula in both places is how a measured artifact and a shipped
artifact silently drift apart.

The field is

    S = normalise( magnetic_family_mean )          # rtp, tmi, mag_anom -- independent acquisitions
        * (0.5 + 0.5 * normalise(gravity))         # iso_grav_anom -- an independent potential field

so a fault-like contact expressed in BOTH potential-field families keeps its full weight while a
magnetic-only (or gravity-only) cluster keeps half.  The element-wise-minimum variant was measured
too and is reported in ``data/evidence/crest_variants.json``; the corroboration rule above scored
higher on the blocked instrument, so it is the one that ships.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
from scipy.ndimage import maximum_filter

MAG_LAYERS = ("rtp", "tmi", "mag_anom")
GRAV_LAYER = "iso_grav_anom"


def normalise(f: np.ndarray, foot: np.ndarray) -> np.ndarray:
    """Scale a density field to [0, 1] inside the footprint (NaN outside)."""
    m = np.nanmax(f[foot])
    if np.isfinite(m) and m > 0:
        return np.where(foot, f / m, np.nan)
    return np.where(foot, 0.0, np.nan)


def crest(f: np.ndarray, foot: np.ndarray) -> np.ndarray:
    """Positive-mass 1-px ridge crests of a continuous density field (plateau-safe).

    ``f > 0`` matters: without it every pixel of a zero plateau is its own local maximum.
    """
    f2 = np.nan_to_num(f, nan=-np.inf)
    mx = maximum_filter(f2, size=3, mode="nearest")
    return (f2 >= mx) & (f > 0) & foot


def crossfamily(fields: dict[str, np.ndarray], foot: np.ndarray) -> np.ndarray:
    """Cross-family corroborated cluster density from per-layer KDE fields."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)      # all-NaN columns outside the footprint
        mag = np.nanmean(np.stack([normalise(fields[k], foot) for k in MAG_LAYERS]), axis=0)
    grav = normalise(fields[GRAV_LAYER], foot)
    s = np.where(foot, np.nan_to_num(mag) * (0.5 + 0.5 * np.nan_to_num(grav)), np.nan)
    return np.where(foot, s / np.nanmax(s[foot]), np.nan)


def robust_norm(f: np.ndarray, foot: np.ndarray) -> np.ndarray:
    """Rank-normalise a density field to [0, 1] inside the footprint.

    Different Euler window scales produce fields with different value distributions (a large window
    averages more, so its densities are smoother and lower-variance).  Comparing raw values across
    scales would let one scale dominate; ranks are scale-free.
    """
    v = np.nan_to_num(f, nan=0.0)[foot]
    r = np.argsort(np.argsort(v)) / max(v.size - 1, 1)
    out = np.full(f.shape, np.nan)
    out[foot] = r
    return out


def multiscale_rankmin(fields_by_scale: dict[int, np.ndarray], foot: np.ndarray) -> np.ndarray:
    """Multi-scale conjunction: the *worst* rank a pixel holds across all window scales.

    A fault-like contact whose Euler solutions are stable across window sizes keeps a high rank at
    every scale, so its rank-min stays high; a solution that only appears at one scale -- the
    classic single-window Euler artefact -- is demoted by the other scales.
    """
    ranks = [robust_norm(f, foot) for f in fields_by_scale.values()]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)      # all-NaN rows outside the footprint
        return np.nanmin(np.stack(ranks), axis=0)


def load_fields(work: str | Path, window: int | None = None,
                needed: tuple[str, ...] = MAG_LAYERS + (GRAV_LAYER,)) -> dict[str, np.ndarray]:
    """Load the per-layer KDE fields written by ``run_euler.py``.

    ``work/`` accumulates one ``.npz`` per (layer, window) as experiments are run, and a naive
    ``glob("euler_*.npz")`` lets the alphabetically-last window silently overwrite the others --
    the failure mode is a shipped artifact built from an unmeasured window.  Select the window
    explicitly and refuse to guess when a layer has more than one candidate file.
    """
    work = Path(work)
    fields: dict[str, np.ndarray] = {}
    files: dict[str, list[str]] = {}
    for npz in sorted(work.glob("euler_*.npz")):
        with np.load(npz, allow_pickle=True) as d:
            if window is not None and int(d["window"]) != int(window):
                continue
            layer = str(d["layer"])
            if layer in needed:
                fields[layer] = d["field"].astype(np.float64)
                files.setdefault(layer, []).append(npz.name)
    duplicated = {k: v for k, v in files.items() if len(v) > 1}
    if duplicated:
        raise SystemExit(f"ambiguous Euler fields for {sorted(duplicated)}: {duplicated} -- "
                         f"pass --window to select one scale")
    missing = set(needed) - set(fields)
    if missing:
        raise SystemExit(f"missing Euler fields for {sorted(missing)} at window={window} -- "
                         f"run scripts/run_euler.py first")
    return fields
