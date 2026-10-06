#!/usr/bin/env python3
"""H8 solver: SI = 0 contact Euler clouds with H8's own gates (H4 stays frozen).

Deliberately a local re-implementation of the window loop so that the frozen
``gemsdoe40.contact_euler.SETTINGS`` used by the published H4 artifact are never
modified.  Only the *math* (``spectral_gradients``, ``solve_contact_windows``) is
imported from the frozen module.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from gemsdoe40.contact_euler import FieldConfig, solve_contact_windows
from gemsdoe40.h8_euler import CLOUD_DTYPE, SETTINGS, FamilyConfig


def solve_window(tx: np.ndarray, ty: np.ndarray, tz: np.ndarray, valid: np.ndarray,
                 family: FamilyConfig, window: int, *, family_id: int,
                 stride: int, cell_m: float = 100.0, batch_size: int = 256):
    """One family/window Euler pass. Returns (cloud, totals)."""
    half = window // 2
    if window < 5 or window % 2 != 1 or stride < 1:
        raise ValueError("odd window >= 5 and positive stride required")
    safe = ndimage.minimum_filter(valid.astype(np.uint8), size=window, mode="constant", cval=0).astype(bool)
    lattice = np.zeros(valid.shape, dtype=bool)
    lattice[half:valid.shape[0] - half:stride, half:valid.shape[1] - half:stride] = True
    rows, cols = np.nonzero(safe & lattice)
    offsets = np.arange(-half, half + 1)
    dr, dc = np.meshgrid(offsets, offsets, indexing="ij")
    x, y = dc.ravel() * cell_m, -dr.ravel() * cell_m
    views = [np.lib.stride_tricks.sliding_window_view(a, (window, window)) for a in (tx, ty, tz)]
    output = []
    totals = {"family": family.name, "window": window, "windows_tested": int(len(rows)),
              "accepted": 0, "rejected_condition": 0, "rejected_residual": 0,
              "rejected_depth": 0, "rejected_horizontal": 0, "rejected_depth_se": 0}
    for start in range(0, len(rows), batch_size):
        rr, cc = rows[start:start + batch_size], cols[start:start + batch_size]
        d = np.stack([v[rr - half, cc - half].reshape(len(rr), -1) for v in views], axis=-1)
        fitted = solve_contact_windows(d, x, y)
        pos = fitted["position"]
        depth = pos[:, 2] - family.continuation_m
        cond_ok = np.isfinite(pos).all(axis=1) & (fitted["condition"] <= SETTINGS["max_design_condition"])
        res_ok = fitted["residual"] <= SETTINGS["max_relative_residual"]
        dep_ok = (depth >= SETTINGS["min_effective_depth_m"]) & (depth <= family.max_depth_m)
        xy_ok = np.hypot(pos[:, 0], pos[:, 1]) <= window * cell_m / 2.0
        se_ok = fitted["se"][:, 2] <= np.maximum(SETTINGS["depth_se_floor_m"],
                                                 SETTINGS["max_depth_se_fraction"] * depth)
        totals["rejected_condition"] += int(np.count_nonzero(~cond_ok))
        totals["rejected_residual"] += int(np.count_nonzero(~res_ok))
        totals["rejected_depth"] += int(np.count_nonzero(~dep_ok))
        totals["rejected_horizontal"] += int(np.count_nonzero(~xy_ok))
        totals["rejected_depth_se"] += int(np.count_nonzero(~se_ok))
        good = cond_ok & res_ok & dep_ok & xy_ok & se_ok
        sol = np.zeros(int(good.sum()), dtype=CLOUD_DTYPE)
        sol["row"] = rr[good] - pos[good, 1] / cell_m
        sol["col"] = cc[good] + pos[good, 0] / cell_m
        sol["depth_m"] = depth[good]
        sol["depth_se_m"] = fitted["se"][good, 2]
        sol["residual"] = fitted["residual"][good]
        sol["condition"] = fitted["condition"][good]
        sol["rank"] = fitted["rank"][good]
        sol["offset"] = fitted["offset"][good]
        sol["window"] = window
        sol["family_id"] = family_id
        output.append(sol)
    cloud = np.concatenate(output) if output else np.empty(0, dtype=CLOUD_DTYPE)
    totals["accepted"] = int(len(cloud))
    totals["rank2_accepted"] = int(np.count_nonzero(cloud["rank"] == 2))
    return cloud, totals


def field_config(family: FamilyConfig, window: int) -> FieldConfig:
    return FieldConfig(family.name, family.continuation_m, family.derivative_order,
                       (window, window), family.edge_guard_m, family.max_depth_m)
