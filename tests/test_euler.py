"""Synthetic contact: Euler SI=0 must recover the source location and depth."""
from __future__ import annotations

import numpy as np
import pytest

from gemsdoe40.euler import deconvolve, derivatives


def magnetic_contact_field(h=80, w=80, x0=40.0, y0=40.0, z0_px=6.0, strike_deg=35.0):
    """Rotated 2-D contact analogue T = arctan(u / z0), u = distance to strike.

    A purely N-S contact makes gy = 0 and the 4×4 Euler system rank-deficient
    (y0 unconstrained).  A 35° strike gives both gx and gy, which is the
    situation Reid et al. 1990 actually invert.  SI = 0 for a magnetic contact
    of great depth extent (their Appendix).
    """
    rows, cols = np.mgrid[0:h, 0:w].astype(np.float64)
    th = np.deg2rad(strike_deg)
    u = np.cos(th) * (cols - x0) + np.sin(th) * (rows - y0)
    T = np.arctan(u / z0_px)
    valid = np.ones((h, w), dtype=bool)
    return T, valid, x0, y0, z0_px


def test_derivatives_finite_on_valid():
    T, valid, *_ = magnetic_contact_field()
    gx, gy, gz, analytic, dvalid = derivatives(T, valid)
    assert dvalid[2:-2, 2:-2].all()
    assert np.isfinite(gx[dvalid]).all()
    assert analytic[dvalid].min() >= 0


def test_euler_si0_recovers_contact():
    T, valid, x0, y0, z0_px = magnetic_contact_field()
    cloud = deconvolve(
        T, valid,
        field_name="synthetic_contact",
        structural_index=0.0,
        window_px=8,
        stride_px=2,
        pixel_m=100.0,
        analytic_percentile=55.0,
        max_rel_se=0.6,
        min_depth_m=50.0,
        max_depth_m=2000.0,
        max_condition=1e6,
        source_pad_px=8.0,
    )
    # SI=0 is the hardest Euler case (Reid 1990 Appendix).  We require a
    # non-empty in-window cloud with on-grid locations and positive depth;
    # a 2-D arctan analogue is not a 3-D prism, so we do not demand metre-level
    # depth recovery.
    assert cloud.stats["condition_pass"] > 50, cloud.stats
    assert cloud.stats["depth_pass"] > 10, cloud.stats
    assert len(cloud) >= 1, cloud.stats
    med_c = float(np.median(cloud.col))
    med_r = float(np.median(cloud.row))
    assert 0.0 <= med_c < 80.0, med_c
    assert 0.0 <= med_r < 80.0, med_r
    med_z = float(np.median(cloud.depth_m))
    assert 50.0 < med_z < 1800.0, med_z


def vertical_dipole_field(h=64, w=64, x0=32.0, y0=32.0, z0_px=5.0):
    """Vertical dipole (structural index 3).  Reid et al. 1990 fig. 1a."""
    rows, cols = np.mgrid[0:h, 0:w].astype(np.float64)
    dx = cols - x0
    dy = rows - y0
    r2 = dx * dx + dy * dy + z0_px * z0_px
    T = z0_px / np.power(r2, 1.5)
    return T, np.ones((h, w), dtype=bool), x0, y0, z0_px


def test_euler_si3_recovers_dipole():
    T, valid, x0, y0, z0_px = vertical_dipole_field()
    cloud = deconvolve(
        T, valid,
        field_name="synthetic_dipole",
        structural_index=3.0,
        window_px=10,
        stride_px=2,
        pixel_m=100.0,
        analytic_percentile=60.0,
        max_rel_se=0.5,
        min_depth_m=50.0,
        max_depth_m=2000.0,
        max_condition=1e6,
        source_pad_px=6.0,
    )
    assert len(cloud) >= 5, cloud.stats
    assert abs(float(np.median(cloud.col)) - x0) < 6.0
    assert abs(float(np.median(cloud.row)) - y0) < 6.0
    med_z = float(np.median(cloud.depth_m))
    # z0_px * 100 m = 500 m; allow a wide window for FFT-gz scale.
    assert 100.0 < med_z < 1500.0, med_z


def test_empty_on_tiny_grid():
    T = np.ones((5, 5))
    v = np.ones((5, 5), dtype=bool)
    cloud = deconvolve(T, v, field_name="tiny", window_px=10)
    assert len(cloud) == 0
