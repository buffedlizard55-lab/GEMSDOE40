"""Preregistered synthetic verification for the H8 ASA/local-wavenumber estimator.

Sources: the contact anomaly is built from Im(log(x-x0 + i(z-z0))) with a deep
regional compensation so the profile is band-limited like real leveled TMI;
the analytic-signal envelope is the Nabighian (1972) form
|C|/sqrt((x-x0)^2 + z0^2); the depth estimate z = sqrt(A/|lambda_min(A)|) is
the contact-calibrated local-wavenumber depth of the SPI family (Thurston &
Smith 1997), implemented with the registered supergaussian spectral taper.
"""
from __future__ import annotations

import numpy as np
import pytest

from gemsdoe40 import asa_spi

PAD = 48
RES = asa_spi.RESOLUTION_M
Z_REG_M = 12000.0


def contact_grid(shape, col0_px, z0_m, c=100.0):
    """Regional-removed 2-D contact: C*[atan2(-z0, x-x0) - atan2(-z_reg, x-x0)]."""
    x = (np.arange(shape[1])[None, :] - col0_px) * RES
    prof = c * (np.arctan2(-float(z0_m), x) - np.arctan2(-Z_REG_M, x))
    return np.broadcast_to(prof, shape).copy()


def pole_grid(shape, row0_px, col0_px, z0_m, m=1.0):
    """phi = m/r, r = sqrt((x-x0)^2+(y-y0)^2+z0^2); analytic signal decays as 1/r^2."""
    x = (np.arange(shape[1])[None, :] - col0_px) * RES
    y = (np.arange(shape[0])[:, None] - row0_px) * RES
    r = np.sqrt(x * x + y * y + float(z0_m) ** 2)
    return m / r


def asa_depth(field):
    return asa_spi.analytic_signal_depth(field, pad_cells=PAD)


@pytest.mark.parametrize("z0_m,tol", [(300.0, 0.50), (800.0, 0.20), (1500.0, 0.20)])
def test_s1_contact_depth_and_location(z0_m, tol):
    shape = (257, 257)
    col0 = 128
    out = asa_depth(contact_grid(shape, col0, z0_m))
    a, z = out["a_smooth"], out["z"]
    row = shape[0] // 2
    peak_col = int(np.argmax(a[row, :]))
    assert abs(peak_col - col0) <= 1
    # Nabighian envelope near the peak (taper attenuates the sharpest ridges)
    x = (np.arange(shape[1]) - col0) * RES
    envelope = 100.0 / np.sqrt(x * x + z0_m * z0_m)
    near = slice(col0 - 6, col0 + 7)
    rel = np.abs(a[row, near] - envelope[near]) / envelope[near]
    assert float(np.median(rel)) < 0.15
    assert a[row, peak_col] / envelope[col0] == pytest.approx(1.0, abs=0.15)
    est = float(z[row, peak_col])
    assert abs(est - z0_m) / z0_m <= tol, f"z0={z0_m}, estimated {est}"


def test_s2_dyke_two_edges():
    shape = (257, 257)
    c0 = 116.0
    half_w = 10.0  # 2 km wide dyke
    z0 = 500.0
    field = contact_grid(shape, c0 - half_w, z0) - contact_grid(shape, c0 + half_w, z0)
    out = asa_depth(field)
    a, z = out["a_smooth"], out["z"]
    row = shape[0] // 2
    prof = a[row, :]
    dil = asa_spi.ndimage.maximum_filter(prof, size=3)
    peaks = np.nonzero((prof == dil) & (prof > 0.3 * prof.max()))[0]
    found = []
    for target in (c0 - half_w, c0 + half_w):
        near = peaks[np.abs(peaks - target) <= 3.0]
        assert near.size >= 1, f"no peak near edge {target}: {peaks}"
        found.append(int(near[np.argmax(prof[near])]))
    for pcol in found:
        est = float(z[row, pcol])
        assert abs(est - z0) / z0 <= 0.50, f"edge {pcol}: est {est} vs {z0}"


def test_s3_si_scaling_pole():
    shape = (257, 257)
    z0 = 1000.0
    out = asa_depth(pole_grid(shape, 128, 128, z0))
    a, z = out["a_smooth"], out["z"]
    peak = np.unravel_index(np.argmax(a), a.shape)
    assert abs(peak[0] - 128) <= 1 and abs(peak[1] - 128) <= 1
    est = float(z[peak])
    # A decays as 1/r^2 (n=1)  =>  z_est = z_true / sqrt(2)
    target = z0 / np.sqrt(2.0)
    assert abs(est - target) / target <= 0.15, f"pole depth est {est} vs {target}"


def test_s4_noise_robust():
    rng = np.random.default_rng(8)
    shape = (257, 257)
    z0 = 800.0
    clean = contact_grid(shape, 128, z0)
    noisy = clean + rng.normal(0.0, 0.02 * np.ptp(clean), size=shape)
    out = asa_depth(noisy)
    a, z = out["a_smooth"], out["z"]
    row = shape[0] // 2
    peak_col = int(np.argmax(a[row, :]))
    assert abs(peak_col - 128) <= 1
    med = float(np.median(z[row, max(0, peak_col - 2):peak_col + 3]))
    assert abs(med - z0) / z0 <= 0.30, f"noisy depth {med}"


def test_s5_exclusion_margin():
    shape = (64, 64)
    field = contact_grid(shape, 32, 500.0)
    valid = np.ones(shape, dtype=bool)
    valid[:, 50:] = False  # invalid region on the right flank
    far = asa_spi.far_from_invalid(valid, margin_m=1500.0)
    assert not far[:, 50:].any()
    assert not far[:, 36:50].any()  # within 1500 m of the invalid edge
    assert far[16:48, 16:35].all()  # interior block farther than 1500 m from both
    out = asa_spi.analytic_signal_depth(field, pad_cells=16)
    row, col, z, stats = asa_spi.select_peaks(
        out["a_smooth"], out["z"], valid, margin_m=1500.0, quantile=0.5,
        depth_bounds_m=(100.0, 5000.0))
    assert row.size > 0
    assert (col <= 35).all()


def test_catalogue_ramp_and_support():
    labels = np.zeros((20, 20), dtype=np.int8)
    labels[10, 10] = 1
    ramp = asa_spi.catalogue_ramp(labels, zero_m=150.0, full_m=450.0)
    assert ramp[10, 10] == 0.0
    assert ramp[10, 11] == 0.0  # 100 m away, inside zero zone
    assert ramp[10, 14] == pytest.approx((400.0 - 150.0) / 300.0, abs=1e-6)
    assert ramp[0, 0] == 1.0


def test_dense_weights_ridge():
    """Dense construction: the contact ridge carries weight; quiet flanks do not."""
    shape = (257, 257)
    z0 = 800.0
    field = contact_grid(shape, 128, z0)
    out = asa_depth(field)
    a, z = out["a_smooth"], out["z"]
    valid = np.ones(shape, dtype=bool)
    depth_ok = np.isfinite(z) & (z >= 100) & (z <= 5000)
    thresh = np.quantile(a[depth_ok], 0.90)
    sol = depth_ok & (a >= thresh)
    assert sol.sum() > 100
    zc = z.copy()  # perfectly persistent twin
    w, persist = asa_spi.dense_family_weights(z, sol, a, zc)
    assert persist.any()
    assert (w > 0).sum() > 50
    ridge_mass = w[:, 120:137].sum()
    flank_mass = w[:, :60].sum() + w[:, 197:].sum()
    assert ridge_mass > flank_mass
    kde = asa_spi.ndimage.gaussian_filter(w, sigma=2.0, truncate=3.0)
    assert kde[:, 120:137].mean() > kde[:, :60].mean()


def test_weights_and_kde_mass():
    row = np.array([10.0, 10.4, 30.0])
    col = np.array([10.0, 10.6, 30.0])
    z = np.array([300.0, 320.0, 2000.0])
    persistent = np.array([True, True, False])
    corroborated = np.array([True, False, False])
    w = asa_spi.cloud_weights(row, col, z, persistent, corroborated)
    assert w[2] == 0.0  # isolated point has no neighbours -> n/(n+3) = 0
    assert w[0] > w[1] > 0.0  # corroborated point outranks its uncorroborated twin
    kde = asa_spi.splat_kde(row, col, w, (40, 40), sigma_px=2.0, truncate=3.0)
    assert kde.shape == (40, 40)
    assert float(kde.sum()) == pytest.approx(float(w.sum()), rel=0.01)
    assert kde.argmax() == 10 * 40 + 10  # mass centred on the strongest cluster
