"""Synthetic verification of the H45 Euler solver.

Every test builds a field whose source position and structural index are known *by
construction*, runs the solver, and asserts the recovery.  These are the checks that
separate a working implementation from a plausible-looking one.

Forward models
--------------
``line_mass_grid`` integrates the exact 2-D line-mass kernel numerically, so no closed-form
gravity/magnetic expression has to be trusted.  The analytic-index models
(``point_mass``, ``line_mass``, ``thin_sheet``) use fields that are homogeneous of a known
degree, which is the definition Reid et al. (1990) use for the structural index.
"""
from __future__ import annotations

import numpy as np
import pytest

from gemsdoe40.euler_h45 import (Family, SolveConfig, SolutionCloud, run_family,
                                 solve_euler_grid, spectral_derivatives)

DX = 100.0


def _grid(n=201, dx=DX):
    """Return (east, north) coordinate grids matching a north-up raster.

    The competition rasters have ``dy = -100``, i.e. row index increases *southward*.
    The synthetic grids use the same convention, so ``north = -(row - n//2) * dx``.
    """
    east = (np.arange(n) - n // 2) * dx
    north = -(np.arange(n) - n // 2) * dx
    return np.meshgrid(east, north)


# ----------------------------------------------------------------------------------
# analytic homogeneous fields
# ----------------------------------------------------------------------------------
def point_mass(x0, y0, z0, c=1.0e6):
    """Vertical gravity of a point mass: homogeneous of degree 3 -> SI(G) = 2."""
    e, n = _grid()
    rho2 = (e - x0) ** 2 + (n - y0) ** 2
    return c * z0 / (rho2 + z0 ** 2) ** 1.5


def line_mass(x0, z0, c=1.0e6):
    """Vertical gravity of an infinite horizontal line: homogeneous of degree 2 -> SI = 1."""
    e, _ = _grid()
    return c * z0 / ((e - x0) ** 2 + z0 ** 2)


def horizontal_line_of_dipoles(x0, z0, c=1.0e12):
    """Vertical field of a 2-D line of dipoles (horizontal cylinder / thin bed fault edge).

    Reid & Thurston (2014) Table 1: "horizontal line (cylinder)" is a *line of dipoles*
    with SI(M) = 2, one more than the line-of-poles (thin sheet edge) value of 1, because
    the magnetic field is one vertical derivative of the gravity analogue.
    """
    e, _ = _grid()
    rho2 = (e - x0) ** 2
    return c * (z0 ** 2 - rho2) / (rho2 + z0 ** 2) ** 2


def slab_step(x0, top, bottom, half_width=40000.0, drho=300.0, g=6.674e-11):
    """Vertical gravity of a 2-D **finite** step, integrated analytically.

    Geometry: a block whose right edge is the vertical plane ``east = x0``, occupying
    ``x0 - half_width <= east <= x0`` and ``top <= depth <= bottom``.  This is the
    "finite contact/fault" of Reid & Thurston (2014) Table 1, whose gravity SI is -1.

    Both integrals of the line-mass kernel are closed form::

        int_{d1}^{d2} z / (w^2 + z^2) dz    = 0.5 ln((w^2 + d2^2) / (w^2 + d1^2))
        int      ln(w^2 + c^2) dw           = w ln(w^2 + c^2) - 2w + 2c arctan(w / c)

    so the model is exact and cheap: no quadrature error leaks into a test whose purpose is
    to validate the solver rather than the forward operator.
    """
    e, _ = _grid()

    def horizontal_log_integral(c):
        """int_a^b ln((east - x')^2 + c^2) dx', with a = x0 - half_width and b = x0."""
        wa = e - (x0 - half_width)
        wb = e - x0

        def primitive(w):
            return w * np.log(w * w + c * c) - 2.0 * w + 2.0 * c * np.arctan2(w, c)

        return primitive(wa) - primitive(wb)

    inner = 0.5 * (horizontal_log_integral(bottom) - horizontal_log_integral(top))
    return 2.0 * g * drho * inner * 1e5   # SI -> mGal


# ----------------------------------------------------------------------------------
# tests
# ----------------------------------------------------------------------------------
def _at(res, iy=100, ix=100):
    """Index of the stride-lattice element whose window centre is grid cell (iy, ix)."""
    jy = int(np.where(res["row_idx"] == iy)[0][0])
    jx = int(np.where(res["col_idx"] == ix)[0][0])
    return jy, jx


def _solve(field, si, window=15, stride=1, kappa=0.25):
    ny, nx = field.shape
    valid = np.ones((ny, nx), bool)
    valid[:8] = False           # exercise the invalid-cell machinery on a clean model
    valid[-8:] = False
    fld, fx, fy, fz = spectral_derivatives(field, valid, DX, pad=48)
    n_valid = np.full((ny, nx), float(window * window))
    cfg = SolveConfig(dx=DX, min_along_strike_power=kappa, max_depth_m=1e6,
                      max_rel_residual=1.0, max_depth_se_frac=1e6,
                      max_abs_offset_frac=1e6, min_window_valid_frac=0.0)
    return solve_euler_grid(fx, fy, fz, fld, si, window, stride, cfg, n_valid)


def test_east_north_down_sign_conventions():
    """d/dx must be the east derivative and d/dy the north derivative on a north-up grid."""
    e, n = _grid(201)
    field = np.exp(-((e - 3000.0) ** 2 + (n + 2000.0) ** 2) / (2 * 3000.0 ** 2))
    valid = np.ones(field.shape, bool)
    _, fx, fy, fz = spectral_derivatives(field, valid, DX, pad=48)
    # east of the peak the field decreases: dF/dx < 0 there
    iy, ix = np.unravel_index(np.argmax(field), field.shape)
    assert fx[iy, min(ix + 40, 200)] < 0.0
    # north of the peak is a *smaller* row index
    assert fy[max(iy - 40, 0), ix] < 0.0
    # a positive anomaly above the plane has a positive downward derivative
    assert fz[iy, ix] > 0.0


def test_point_mass_recovers_depth_and_position_at_si2():
    z0, x0, y0 = 1500.0, -400.0, 700.0
    field = point_mass(x0, y0, z0)
    res = _solve(field, si=2.0, window=15, stride=1)
    cy, cx = _at(res)
    got_z = res["z0"][cy, cx]
    assert abs(got_z - z0) / z0 < 0.05, f"depth {got_z} vs {z0}"
    east = res["centre_col"][cy, cx] * DX - 100 * DX + res["x0"][cy, cx]
    north = -(res["centre_row"][cy, cx] * DX - 100 * DX) + res["y0"][cy, cx]
    assert abs(east - x0) < 250.0, f"east {east} vs {x0}"
    assert abs(north - y0) < 250.0, f"north {north} vs {y0}"


def test_line_mass_recovers_depth_at_si1():
    z0, x0 = 1200.0, 600.0
    field = line_mass(x0, z0)
    res = _solve(field, si=1.0, window=15, stride=1)
    cy, cx = _at(res)
    assert abs(res["z0"][cy, cx] - z0) / z0 < 0.05
    assert abs(res["x0"][cy, cx] - x0) < 250.0
    # the along-strike coordinate is unidentifiable, so the regulariser must pin it near 0
    assert abs(res["y0"][cy, cx]) < 150.0


def test_horizontal_line_of_dipoles_recovers_depth_at_si2():
    """A line of dipoles is one derivative above a line of poles: SI(M) = 2, not 1."""
    z0, x0 = 900.0, -300.0
    field = horizontal_line_of_dipoles(x0, z0)
    res = _solve(field, si=2.0, window=13, stride=1)
    cy, cx = _at(res)
    assert abs(res["z0"][cy, cx] - z0) / z0 < 0.08
    assert abs(res["x0"][cy, cx] - x0) < 250.0
    # solving the same model at the wrong index must move the depth estimate substantially
    wrong = _solve(field, si=1.0, window=13, stride=1)
    assert abs(wrong["z0"][cy, cx] - res["z0"][cy, cx]) > 0.10 * z0


def test_straight_contact_is_pinned_along_strike_not_across():
    """For a y-invariant field the along-strike offset must shrink to ~0 while the
    across-strike offset must still be recovered."""
    z0, x0 = 1000.0, 800.0
    field = line_mass(x0, z0)
    res = _solve(field, si=1.0, window=15, stride=1)
    cy, cx = _at(res)
    assert abs(res["y0"][cy, cx]) < 150.0          # along strike: pinned
    assert abs(res["x0"][cy, cx] - x0) < 250.0     # across strike: free


def test_gravity_finite_step_edge_is_recovered_in_position():
    """Reid & Thurston (2014): the gravity SI of a finite step is -1.  The reported depth is
    biased for a finite step (they warn it needs a generalised formulation), but the *edge
    position* must still be recovered.  This asserts position, not depth."""
    top, bottom, x0 = 200.0, 2200.0, 0.0
    field = slab_step(x0, top, bottom)
    res = _solve(field, si=-1.0, window=25, stride=1)
    cy, cx = _at(res)
    east = res["centre_col"][cy, cx] * DX - 100 * DX + res["x0"][cy, cx]
    assert abs(east - x0) < 400.0, f"edge at {east}, expected {x0}"


def test_si_zero_offset_form_matches_reid_equation_2():
    """At SI=0 the solver must use Reid et al. eq. (2) with the free offset A, i.e. the
    fourth design column is 1, not 0.  Checked by construction and by a scaled model:
    multiplying the field by a constant must not move the solution."""
    z0, x0 = 1000.0, -500.0
    f1 = line_mass(x0, z0)
    r1 = _solve(f1, si=0.0, window=15, stride=1)
    r2 = _solve(f1 * 3.7 + 1234.0, si=0.0, window=15, stride=1)
    cy, cx = _at(r1)
    assert abs(r1["x0"][cy, cx] - r2["x0"][cy, cx]) < 1.0
    assert abs(r1["z0"][cy, cx] - r2["z0"][cy, cx]) < 1.0


def test_relative_residual_is_scale_free():
    field = line_mass(0.0, 1000.0)
    r1 = _solve(field, si=1.0, window=15, stride=1)
    r2 = _solve(field * 5.0, si=1.0, window=15, stride=1)
    # compared with a loose absolute floor: in windows with essentially no signal the
    # misfit is ~1e-6 and the comparison is dominated by floating-point round-off
    np.testing.assert_allclose(r1["rel_resid"], r2["rel_resid"], rtol=1e-4, atol=1e-6)


def test_run_family_quality_gates_reduce_the_cloud():
    field = point_mass(0.0, 0.0, 1800.0)
    valid = np.ones(field.shape, bool)
    fam = Family(name="synth", kind="gravity", si=2.0, window=15, stride=2)
    cfg = SolveConfig(dx=DX, max_depth_m=4000.0)
    cloud = run_family(field, valid, fam, cfg)
    assert isinstance(cloud, SolutionCloud)
    assert len(cloud) > 0
    assert cloud.meta["counts"]["windows"] > len(cloud)
    assert np.all(cloud.depth >= cfg.min_depth_m)
    assert np.all(cloud.depth <= cfg.max_depth_m)
    assert np.all(np.isfinite(cloud.depth_se))


def test_spectral_derivative_of_a_constant_is_zero():
    field = np.full((121, 121), 7.5)
    valid = np.ones(field.shape, bool)
    _, fx, fy, fz = spectral_derivatives(field, valid, DX, pad=32)
    assert np.abs(fx).max() < 1e-9
    assert np.abs(fy).max() < 1e-9
    assert np.abs(fz).max() < 1e-9
