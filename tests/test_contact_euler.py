"""Analytic tests independent of prior Euler code and the competition proxy."""
from dataclasses import replace

import numpy as np
import pytest

from gemsdoe40.contact_euler import (
    CLOUD_DTYPE, FIELDS, bilinear_splat, cluster_weights, combine_families,
    family_kde, solution_cloud, solve_contact_windows, spectral_gradients,
)


def contact(x, y, normal, location=125.0, depth=900.0, angular=2.0, logarithmic=0.7):
    """Independent harmonic F=C atan(u/q)+L log(sqrt(u²+q²)), q=depth-z.

    On z=0, u*Fn-depth*Fz=L. This exercises Reid's nonzero offset,
    rather than testing a zero-offset simplified source against itself.
    """
    nx, ny = normal
    u = nx * x + ny * y - location
    r2 = u*u + depth*depth
    fn = (angular * depth + logarithmic * u) / r2
    fz = (angular * u - logarithmic * depth) / r2
    return np.stack((nx * fn, ny * fn, fz), axis=-1)


def coordinates(size=25, cell=100):
    offsets = (np.arange(size) - size//2) * cell
    x, row = np.meshgrid(offsets, offsets)
    return x, -row


@pytest.mark.parametrize("angle", [0, 15, 45, 90, 130, 175])
@pytest.mark.parametrize("amplitude", [1e-6, 1, 1e5])
def test_rotated_contact_recovers_nonzero_offset_depth_and_crossstrike(angle, amplitude):
    x, y = coordinates()
    normal = np.array([np.cos(np.deg2rad(angle)), np.sin(np.deg2rad(angle))])
    d = contact(x, y, normal) * amplitude
    fit = solve_contact_windows(d.reshape(1, -1, 3), x, y)
    assert fit["rank"][0] == 2
    np.testing.assert_allclose(fit["position"][0, :2] @ normal, 125, atol=1e-7)
    np.testing.assert_allclose(fit["position"][0, 2], 900, atol=1e-7)
    np.testing.assert_allclose(fit["offset"][0], 0.7 * amplitude, rtol=1e-9)
    assert fit["residual"][0] < 1e-10
    assert fit["se"][0, 2] < 1e-6


def test_nonzero_offset_is_not_optional_for_general_contact():
    x, y = coordinates()
    d = contact(x, y, (1, 0), logarithmic=1.8).reshape(-1, 3)
    rhs = x.ravel()*d[:, 0] + y.ravel()*d[:, 1]
    wrong = np.linalg.lstsq(d[:, [0, 2]], rhs, rcond=None)[0]
    correct = solve_contact_windows(d[None], x, y)
    assert abs(wrong[1] - 900) > 100
    assert correct["position"][0, 2] == pytest.approx(900)


def test_three_dimensional_identifiable_contact_intersection():
    x, y = coordinates()
    # Shallow intersecting contacts have two well-separated horizontal modes.
    # Deep, broad contacts can legitimately be reduced by the frozen 0.10 rule.
    d = (contact(x, y, (1, 0), location=70, depth=400, logarithmic=0.3)
         + contact(x, y, (0, 1), location=-80, depth=400, angular=1.5, logarithmic=-0.5))
    fit = solve_contact_windows(d.reshape(1, -1, 3), x, y)
    assert fit["rank"][0] == 3
    np.testing.assert_allclose(fit["position"][0], [70, -80, 400], atol=1e-7)
    assert fit["offset"][0] == pytest.approx(-0.2)


def test_translation_changes_only_identifiable_horizontal_location():
    x, y = coordinates()
    n = np.array([0.6, 0.8])
    d = contact(x, y, n).reshape(1, -1, 3)
    fit = solve_contact_windows(d, x + 125000, y - 87000)
    assert fit["position"][0, 2] == pytest.approx(900, abs=1e-6)
    assert fit["position"][0, :2] @ n == pytest.approx(125 + 125000*.6 - 87000*.8)


def test_constant_fields_do_not_invent_a_depth():
    x, y = coordinates()
    d = np.ones((2, x.size, 3))
    fit = solve_contact_windows(d, x, y)
    assert np.isnan(fit["position"]).all()
    assert np.isinf(fit["condition"]).all()
    with pytest.raises(ValueError, match="non-finite"):
        solve_contact_windows(d * np.nan, x, y)


def test_observation_height_removed_once_not_added():
    x, y = coordinates(61)
    d = contact(x, y, (1, 0), location=0, depth=1100)
    conf = replace(FIELDS[0], continuation_m=200)
    cloud, info = solution_cloud(d[:, :, 0], d[:, :, 1], d[:, :, 2], np.ones(x.shape, bool), conf, 15)
    assert info["accepted"] > 0
    np.testing.assert_allclose(cloud["depth_m"], 900, atol=0.01)


def test_fourier_north_and_downward_signs():
    n = 257
    r, c = np.mgrid[:n, :n]
    k = 2*np.pi/64
    f = np.cos(k*c) + 0.6*np.cos(k*r)
    tx, ty, tz = spectral_gradients(f, np.ones(f.shape, bool), pad=128)
    center = np.s_[64:193, 64:193]
    for got, expect in ((tx, -k/100*np.sin(k*c)), (ty, 0.6*k/100*np.sin(k*r)), (tz, k/100*f)):
        assert np.corrcoef(got[center].ravel(), expect[center].ravel())[0, 1] > 0.999
        assert np.std(got[center])/np.std(expect[center]) == pytest.approx(1, rel=0.02)


def cloud_fixture(depth=600):
    cloud = np.zeros(6, dtype=CLOUD_DTYPE)
    cloud["row"] = [20, 21, 20, 19, 21, 19]
    cloud["col"] = [20, 20, 21, 20, 21, 19]
    cloud["depth_m"] = depth
    cloud["depth_se_m"] = depth * .05
    cloud["residual"] = .1
    cloud["window"] = [15, 25, 15, 25, 15, 25]
    return cloud


def test_consistent_shallow_crossscale_clusters_receive_more_weight():
    shallow = cluster_weights(cloud_fixture())
    deep = cluster_weights(cloud_fixture(3000))
    assert np.all(shallow["weight"] > deep["weight"])
    assert np.all(shallow["weight"] > 0)
    bad = cloud_fixture()
    bad["depth_m"] = [100, 2000, 4000, 8000, 16000, 32000]
    assert not cluster_weights(bad)["weight"].any()
    one_scale = cloud_fixture(); one_scale["window"] = 15
    assert not cluster_weights(one_scale)["weight"].any()


def test_bilinear_splat_conserves_mass_without_row_col_swap():
    field = bilinear_splat(np.array([2.25]), np.array([5.75]), np.array([8]), (10, 10))
    assert field.sum() == pytest.approx(8)
    assert field[2, 5] == pytest.approx(1.5)
    assert field[2, 6] == pytest.approx(4.5)
    assert field[3, 5] == pytest.approx(.5)
    assert field[3, 6] == pytest.approx(1.5)


def test_continuous_kde_has_finite_range_and_exact_mask():
    foot = np.ones((41, 41), bool); foot[0] = False
    known = np.zeros_like(foot); known[20, 20] = True
    cloud = cluster_weights(cloud_fixture())
    kde, receipt = family_kde(cloud, foot)
    p = combine_families(kde, kde, foot, known)
    assert receipt["clustered_solutions"] == 6
    assert np.isfinite(p).all() and p.min() == 0 and p.max() == 1
    assert np.unique(p).size > 20
    assert p[20, 20] == 0 and not p[0].any()
    # Exact known-cell mask, not a >=100 m buffer.
    assert p[20, 21] > 0
    with pytest.raises(ValueError, match="no clustered"):
        family_kde(np.empty(0, dtype=CLOUD_DTYPE), foot)


def test_missing_data_window_not_used():
    x, y = coordinates(51)
    d = contact(x, y, (1, 0), depth=1200)
    valid = np.ones(x.shape, bool); valid[25, 25] = False
    cloud, info = solution_cloud(d[:, :, 0], d[:, :, 1], d[:, :, 2], valid, FIELDS[0], 15)
    assert info["windows_tested"] < 64
    assert len(cloud) <= info["windows_tested"]
