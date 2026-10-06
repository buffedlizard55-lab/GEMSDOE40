import numpy as np
import pytest
import rasterio

from gemsdoe40.euler_h4 import (
    EulerCloud,
    PairCloud,
    VoxelCloud,
    mutual_nearest_pairs,
    project_pair_cloud_to_kde,
    solve_euler_cloud,
    voxelize_and_weight,
)


CELL = 100.0
H = W = 40
TRANSFORM = rasterio.Affine(CELL, 0.0, 0.0, 0.0, -CELL, H * CELL)


def _synthetic_homogeneous_field(si: float):
    rows, cols = np.mgrid[0:H, 0:W]
    x = (cols + 0.5) * CELL
    y = H * CELL - (rows + 0.5) * CELL
    x0, y0, z0 = 2030.0, 2000.0, 700.0
    u, v, w = x - x0, y - y0, -z0
    radius = np.sqrt(u * u + v * v + w * w)
    if si == 0.0:
        # A bounded degree-zero homogeneous field with all three spatial directions observable.
        T = u / radius + 4.0
        tx = (v * v + w * w) / radius**3
        ty = -u * v / radius**3
        tz = -u * w / radius**3
    elif si == -1.0:
        # Degree +1 field (SI=-1) plus a regional constant B=5.
        T = radius + 5.0
        tx = u / radius
        ty = v / radius
        tz = w / radius
    else:
        raise AssertionError(si)
    return T.astype(np.float32), tx.astype(np.float32), ty.astype(np.float32), tz.astype(np.float32), x0, y0, z0


@pytest.mark.parametrize("si", [0.0, -1.0])
def test_local_euler_solver_recovers_synthetic_3d_homogeneous_source(si):
    T, tx, ty, tz, x0, y0, z0 = _synthetic_homogeneous_field(si)
    cloud, counts = solve_euler_cloud(
        T,
        tx,
        ty,
        tz,
        np.ones_like(T, dtype=bool),
        transform=TRANSFORM,
        si=si,
        field=f"synthetic-si-{si:g}",
        window_px=10,
        stride_px=2,
        chunk_window_rows=4,
    )
    expected_row = (y0 - TRANSFORM.f) / TRANSFORM.e - 0.5
    expected_col = (x0 - TRANSFORM.c) / TRANSFORM.a - 0.5
    distance = np.hypot(cloud.row - expected_row, cloud.col - expected_col)
    near = distance < 1.0
    assert counts["windows_conditioned"] > 0
    assert near.sum() > 2
    assert np.median(cloud.depth_m[near]) == pytest.approx(z0, abs=1.0)
    assert np.median(cloud.relative_residual[near]) < 1e-5
    assert np.median(distance[near]) < 0.05


def _cloud(rows, cols, depths, residual=None, offset=None, si=0.0, field="test"):
    n = len(rows)
    return EulerCloud(
        row=np.asarray(rows, dtype=np.float32),
        col=np.asarray(cols, dtype=np.float32),
        depth_m=np.asarray(depths, dtype=np.float32),
        relative_residual=np.asarray(residual if residual is not None else [0.02] * n, dtype=np.float32),
        condition=np.full(n, 20.0, dtype=np.float32),
        offset_px=np.asarray(offset if offset is not None else [0.2] * n, dtype=np.float32),
        si=si,
        field=field,
    )


def test_voxelization_removes_overlap_votes_and_requires_tight_3d_neighbours():
    cloud = _cloud(
        rows=[5.0, 5.0, 7.1, 18.0],
        cols=[5.0, 7.1, 6.0, 18.0],
        depths=[500.0, 530.0, 560.0, 500.0],
    )
    footprint = np.ones((24, 24), dtype=bool)
    transform = rasterio.Affine(CELL, 0.0, 0.0, 0.0, -CELL, 2400.0)
    voxels, stats = voxelize_and_weight(cloud, footprint, transform)
    assert stats["in_footprint_solutions"] == 4
    assert stats["unique_source_voxels"] == 4
    assert stats["cluster_supported_voxels"] == 3
    assert np.count_nonzero(voxels.weight > 0) == 3
    assert voxels.neighbour_count[:3].min() == 3
    assert voxels.weight[-1] == 0.0


def test_voxel_base_weight_accepts_bounded_independent_context_factor():
    cloud = _cloud(
        rows=[5.0, 5.0, 7.1],
        cols=[5.0, 7.1, 6.0],
        depths=[500.0, 530.0, 560.0],
    )
    footprint = np.ones((24, 24), dtype=bool)
    transform = rasterio.Affine(CELL, 0.0, 0.0, 0.0, -CELL, 2400.0)
    baseline, _ = voxelize_and_weight(cloud, footprint, transform)
    factor = np.array([0.5, 1.0, 0.75], dtype=np.float32)
    weighted, _ = voxelize_and_weight(cloud, footprint, transform, source_context_factor=factor)
    assert np.allclose(weighted.base_weight / baseline.base_weight, [0.5, 0.75, 1.0])
    with pytest.raises(ValueError):
        voxelize_and_weight(cloud, footprint, transform, source_context_factor=np.array([0.5, 0.75, 1.1]))


def _voxel_cloud(x, y, depth, weight, row=None, col=None):
    x = np.asarray(x, dtype=np.float32)
    y = np.asarray(y, dtype=np.float32)
    depth = np.asarray(depth, dtype=np.float32)
    weight = np.asarray(weight, dtype=np.float32)
    n = len(x)
    z = np.zeros(n, dtype=np.float32)
    return VoxelCloud(
        row=np.asarray(row if row is not None else z, dtype=np.float32),
        col=np.asarray(col if col is not None else z, dtype=np.float32),
        x_m=x,
        y_m=y,
        depth_m=depth,
        base_weight=weight,
        cluster_weight=np.ones(n, dtype=np.float32),
        weight=weight,
        neighbour_count=np.full(n, 3, dtype=np.int16),
        depth_mad_m=np.zeros(n, dtype=np.float32),
    )


def test_mutual_nearest_matching_is_one_to_one_and_applies_3d_gates():
    magnetic = _voxel_cloud([0, 1000, 5000], [0, 0, 0], [500, 800, 900], [1, 0.8, 1], row=[3, 4, 5], col=[3, 10, 20])
    gravity = _voxel_cloud([40, 1100, 5500], [0, 0, 0], [520, 780, 1500], [0.9, 0.7, 1], row=[3, 4, 5], col=[3, 10, 20])
    pairs = mutual_nearest_pairs(magnetic, gravity)
    assert pairs.row.size == 2
    assert pairs.xy_distance_m.max() <= 300
    assert pairs.depth_delta_m.max() <= 300
    assert np.all(pairs.weight > 0)


def test_continuous_kde_zeros_only_exact_known_label_and_nans_outside():
    footprint = np.ones((15, 15), dtype=bool)
    footprint[0, 0] = False
    known = np.zeros((15, 15), dtype=np.int8)
    known[7, 7] = 1
    pairs = PairCloud(
        row=np.array([7.2], dtype=np.float32),
        col=np.array([7.4], dtype=np.float32),
        depth_m=np.array([600.0], dtype=np.float32),
        weight=np.array([1.0], dtype=np.float32),
        xy_distance_m=np.array([50.0], dtype=np.float32),
        depth_delta_m=np.array([10.0], dtype=np.float32),
    )
    prediction, info = project_pair_cloud_to_kde(pairs, footprint, known)
    assert np.isnan(prediction[0, 0])
    assert prediction[7, 7] == 0.0
    assert prediction[7, 8] > 0.0
    assert prediction[footprint].min() >= 0.0
    assert prediction[footprint].max() == pytest.approx(1.0)
    assert info["positive_in_footprint"] > 0


def test_kde_refuses_to_invent_support_when_pair_cloud_is_empty():
    footprint = np.ones((3, 3), dtype=bool)
    labels = np.zeros((3, 3), dtype=np.int8)
    empty = PairCloud(*(np.empty(0, dtype=np.float32) for _ in range(6)))
    with pytest.raises(RuntimeError, match="no positive 3-D pair support"):
        project_pair_cloud_to_kde(empty, footprint, labels)
