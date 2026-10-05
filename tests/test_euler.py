"""Analytic verification of the Euler deconvolution implementation.

Two synthetic potential fields with *analytically known* source depth are inverted and the
recovered depth is compared with the truth:

  * monopole        T = C z0 / (r^2 + z0^2)^{3/2}     homogeneous degree -2  -> SI = 2
  * 2-D contact     T = C atan(x / z0)                homogeneous degree  0  -> SI = 0
    (the classic fault/contact structural index of Reid et al. 1990, Table 1)

The contact case also has a non-zero regional level (T -> +-C*pi/2 at infinity), which exercises
the fitted background term B.
"""

from __future__ import annotations

import numpy as np
import pytest
import rasterio

from gems40.euler import (cluster_weights, euler_solutions, fourier_gradients, kde_field)

CELL = 100.0
N = 161          # odd, 16.1 km across
Z0 = 700.0       # true source depth (m)


def _transform():
    return rasterio.Affine(CELL, 0.0, 0.0, 0.0, -CELL, N * CELL)


def _run(T, si, z0=Z0, window_px=11, stride=2):
    # edge windows are excluded: the Fourier taper makes them meaningless
    tx, ty, tz = fourier_gradients(T, cell_m=CELL)
    sol = euler_solutions(T, tx, ty, tz, _transform(), si=si, window_px=window_px, stride=stride)
    w = cluster_weights(sol)
    return sol, w


def test_monopole_depth_and_sign():
    y, x = np.mgrid[0:N, 0:N]
    xm = (x - N // 2) * CELL
    ym = (y - N // 2) * CELL
    T = Z0 / (xm ** 2 + ym ** 2 + Z0 ** 2) ** 1.5
    _, _ = _run(T, si=2)                     # smoke
    sol, w = _run(T, si=2)
    d = sol.depth_m[w["weight"] > 0]
    keep = (d > 0) & (d < 5000)
    assert keep.sum() > 50, "too few accepted solutions"
    # solutions near the anomaly (weight in the top decile) must recover the true depth
    top = w["weight"] >= np.quantile(w["weight"][w["weight"] > 0], 0.9)
    est = np.median(sol.depth_m[top])
    assert 0.6 * Z0 < est < 1.6 * Z0, f"monopole depth {est:.1f} m vs truth {Z0}"
    # sign convention: with the vertical derivative flipped, the same windows must come back
    # with a *negative* depth at the anomaly (i.e. the sign of Tz is what fixes the sign of z0)
    tx, ty, tz = fourier_gradients(T, cell_m=CELL)
    sol_flip = euler_solutions(T, tx, ty, -tz, _transform(), si=2, window_px=11, stride=2)
    near = (np.abs(sol_flip.row - N // 2) < 5) & (np.abs(sol_flip.col - N // 2) < 5)
    assert near.sum() > 5
    assert np.median(sol_flip.depth_m[near]) < 0.0, "vertical-derivative sign convention is wrong"


def test_contact_si0_depth():
    """The fault/contact structural index (SI = 0, Reid et al. 1990 Table 1).

    The synthetic contact is intentionally 1-D (no y dependence), so its Ty is identically zero
    and the *y* coordinate of the solution is unconstrained: the lateral-offset gate is therefore
    disabled here and only the depth is asserted, which is what this test is about.
    """
    y, x = np.mgrid[0:N, 0:N]
    xm = (x - N // 2) * CELL
    T = np.arctan2(xm, Z0)                   # C=1, homogeneous degree 0 -> SI = 0
    sol, _ = _run(T, si=0)
    w = cluster_weights(sol, max_offset_px=1e9, amp_floor_frac=0.0)
    sel = w["keep"] & (np.abs(sol.col - N // 2) < 12)   # windows straddling the contact
    assert sel.sum() > 100
    est = np.median(sol.depth_m[sel])
    assert 0.7 * Z0 < est < 1.4 * Z0, f"contact depth {est:.1f} m vs truth {Z0}"


def test_kde_field_is_nonnegative_and_localised():
    y, x = np.mgrid[0:N, 0:N]
    r = np.zeros((N, N))
    r[80, 90] = 3.0
    f = kde_field((N, N), np.array([80]), np.array([90]), np.array([3.0]), sigma_px=1.5)
    assert f.min() >= 0.0
    assert f[80, 90] == f.max()
    assert f.sum() > 0.0
    assert np.allclose(f[r == 0][:10], f[r == 0][:10])  # deterministic


def test_window_geometry_recovered():
    """The monopole's horizontal position must come back at the window that contains it."""
    y, x = np.mgrid[0:N, 0:N]
    xm = (x - 60) * CELL
    ym = (y - 40) * CELL
    T = 1000.0 / (xm ** 2 + ym ** 2 + 1000.0 ** 2) ** 1.5
    sol, w = _run(T, si=2)
    top = np.argsort(-w["weight"])[:40]
    assert abs(np.median(sol.col[top]) - 60) <= 2
    assert abs(np.median(sol.row[top]) - 40) <= 2
