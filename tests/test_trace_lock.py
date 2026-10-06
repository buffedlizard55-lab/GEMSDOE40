"""Unit tests for the H8 trace-locked Euler emission (frozen constants only)."""
from __future__ import annotations

import csv
import gzip

import numpy as np
import pytest

from gemsdoe40 import trace_lock as tl
from gemsdoe40.trace_lock import (
    BIN_DIRECTIONS, CLOUD_DTYPE, CLOUD_SHA256, N_RECORDS,
    _shift, build_trace_locked_field, corridor_kernels, load_cloud_csv,
    orientation_bins, ridge_crests, strike_consensus,
)

TRANSFORM = (100.0, 0.0, 243350.0, 0.0, -100.0, 4508550.0)


def make_cloud(rows, cols, depths, weights, families):
    n = len(rows)
    cloud = np.zeros(n, dtype=CLOUD_DTYPE)
    cloud["family"] = families
    cloud["row"] = rows
    cloud["col"] = cols
    cloud["depth_m"] = depths
    cloud["weight"] = weights
    return cloud


def test_pinned_cloud_identity_constants():
    assert CLOUD_SHA256 == "6bed30b224981ef1064256d696c585b76f438011bbae6c0a16c00e6e6f5922ba"
    assert N_RECORDS == 46_656


def test_load_cloud_csv_roundtrip(tmp_path):
    header = ["family", "window_cells", "easting_m", "northing_m", "effective_depth_m",
              "conditional_depth_se_m", "relative_residual", "design_condition", "contact_offset",
              "identifiable_rank", "cluster_weight", "matching_neighbours", "depth_mad_m"]
    rows = []
    # cell (row=10.5, col=20.25): east = 243350 + 100*(col+0.5), north = 4508550 - 100*(row+0.5)
    east = 243350 + 100 * (20.25 + 0.5)
    north = 4508550 - 100 * (10.5 + 0.5)
    rows.append(["tmi", 15, f"{east:.4f}", f"{north:.4f}", "512.5", "40", "0.1", "50", "0", "2",
                 "0.7", "4", "60"])
    rows.append(["iso_grav_anom", 31, f"{east:.4f}", f"{north:.4f}", "900", "90", "0.2", "60", "0",
                 "2", "0.0", "0", "0"])  # zero weight -> dropped
    path = tmp_path / "cloud.csv.gz"
    with gzip.open(path, "wb") as raw:
        import io
        with io.TextIOWrapper(raw, encoding="utf-8", newline="") as text:
            writer = csv.writer(text)
            writer.writerow(header)
            writer.writerows(rows)
    # record count pin enforced only for the real cloud; patch the constant check out
    import gemsdoe40.trace_lock as mod
    original = mod.N_RECORDS
    mod.N_RECORDS = 2
    try:
        cloud = load_cloud_csv(path, TRANSFORM)
    finally:
        mod.N_RECORDS = original
    assert len(cloud) == 1
    assert cloud["family"][0] == 0
    assert cloud["row"][0] == pytest.approx(10.5)
    assert cloud["col"][0] == pytest.approx(20.25)
    assert cloud["depth_m"][0] == pytest.approx(512.5)
    assert cloud["weight"][0] == pytest.approx(0.7)


def test_shift_no_wrap():
    arr = np.zeros((4, 4), dtype=np.float32)
    arr[2, 0] = 1.0
    # out[r, c] = arr[r + dr, c + dc]
    assert _shift(arr, 1, 0)[1, 0] == 1.0
    assert _shift(arr, -1, 0)[3, 0] == 1.0
    corner = np.zeros((4, 4), dtype=np.float32)
    corner[0, 0] = 1.0
    # out[r, c] = arr[r + dr, c + dc]: a +1 shift needs source row -1 -> zero fill, no wrap
    assert _shift(corner, 1, 0).sum() == 0.0
    assert _shift(corner, 0, 1).sum() == 0.0


def test_plateau_thinning_one_pixel():
    consensus = np.zeros((11, 21), dtype=np.float32)
    consensus[4:7, 3:18] = 1.0  # 3-px wide horizontal plateau: gradient along rows
    bins = np.full(consensus.shape, 2, dtype=np.int8)  # gradient along +row
    ridge = ridge_crests(consensus, bins)
    per_column = ridge[:, 10].sum()
    assert per_column == 1  # plateau thinned to a single crest row
    assert ridge[4:7, 10].sum() == 1


def test_orientation_bins_horizontal_stripes():
    field = np.zeros((40, 40), dtype=np.float32)
    field[18:22, :] = 1.0  # horizontal band -> gradient along rows (bin 2)
    bins, aniso = orientation_bins(field, np.ones_like(field, dtype=bool))
    core = bins[20, 5:35]
    assert (core == 2).mean() > 0.9
    assert aniso[20, 20] > tl.ISOTROPY_MIN


def test_corridor_kernels_geometry():
    kernels = corridor_kernels()
    assert set(kernels) == set(BIN_DIRECTIONS)
    k0 = kernels[0]  # gradient along columns -> strike along rows: 17 tall, 7 wide ones
    assert k0.shape == (17, 17)
    assert k0.sum() == 17 * 7
    # across-strike half-width 3 -> columns 5..11 inclusive at the centre row
    assert k0[8, 5] == 1 and k0[8, 11] == 1 and k0[8, 4] == 0 and k0[8, 12] == 0
    assert k0[0, 8] == 1 and k0[16, 8] == 1
    k2 = kernels[2]  # gradient along rows -> strike along columns: 7 tall, 17 wide
    assert k2.sum() == 17 * 7
    assert k2[5, 8] == 1 and k2[11, 8] == 1 and k2[4, 8] == 0 and k2[12, 8] == 0


def test_strike_consensus_accepts_lineament_rejects_scatter():
    rng = np.random.default_rng(7)
    shape = (61, 91)
    footprint = np.ones(shape, dtype=bool)
    # coherent shallow lineament at row 30, both families
    cols = np.arange(10, 70)
    rows_line = np.full(cols.size, 30.0) + rng.normal(0, 0.2, cols.size)
    depths_line = rng.normal(500, 20, cols.size)
    line = make_cloud(np.concatenate([rows_line, rows_line]),
                      np.concatenate([cols, cols]).astype(float),
                      np.concatenate([depths_line, depths_line]),
                      np.ones(2 * cols.size), [0] * cols.size + [1] * cols.size)
    # scattered incoherent cloud; low cluster weights, as the frozen H4 weighting gives
    # deep/scattered solutions in the real cloud
    n = 200
    scat = make_cloud(rng.uniform(0, 60, n), rng.uniform(0, 90, n),
                      rng.uniform(100, 5000, n), np.full(n, 0.05),
                      rng.integers(0, 2, n))
    cloud = np.concatenate([line, scat])
    m, g, _ = tl.family_kdes(cloud, footprint)
    consensus = tl.consensus_field(m, g, footprint, np.zeros(shape, dtype=bool))
    bins, _ = orientation_bins(consensus, footprint)
    ridge = ridge_crests(consensus, bins)
    passes, gamma, stats = strike_consensus(cloud, bins, ridge)
    assert ridge[30, 15:65].mean() > 0.5  # crest detected along the lineament
    assert passes[30, 15:65].mean() > 0.5  # depth consensus keeps it
    assert (gamma[passes] > 0).all() and (gamma[passes] <= 1).all()
    # scattered-only crests must mostly fail the sd gate
    assert (scat["family"] == 0).any() and (scat["family"] == 1).any()
    ms, gs, _ = tl.family_kdes(scat, footprint)
    consensus_s = tl.consensus_field(ms, gs, footprint, np.zeros(shape, dtype=bool))
    bins_s, _ = orientation_bins(consensus_s, footprint)
    ridge_s = ridge_crests(consensus_s, bins_s)
    passes_s, _, stats_s = strike_consensus(scat, bins_s, ridge_s)
    assert stats_s["corridor_passing_cells"] <= 0.05 * max(1, stats_s["ridge_cells"]) + 5


def test_end_to_end_localizes_lineament():
    rng = np.random.default_rng(11)
    shape = (61, 91)
    footprint = np.ones(shape, dtype=bool)
    known = np.zeros(shape, dtype=bool)
    known[30, 40] = True
    cols = np.arange(8, 80)
    rows_line = np.full(cols.size, 30.0) + rng.normal(0, 0.2, cols.size)
    depths_line = rng.normal(400, 15, cols.size)
    cloud = make_cloud(np.concatenate([rows_line, rows_line]),
                       np.concatenate([cols, cols]).astype(float),
                       np.concatenate([depths_line, depths_line]),
                       np.ones(2 * cols.size), [0] * cols.size + [1] * cols.size)
    field, stats = build_trace_locked_field(cloud, footprint, known)
    assert field.shape == shape and field.dtype == np.float32
    assert (field >= 0).all() and (field <= 1).all()
    # normalization happens before exact-known zeroing; if the argmax was a known pixel the
    # stored maximum is the second-highest cell, still bounded by 1
    assert float(field.max()) <= 1.0 and float(field.max()) > 0.9
    assert field[known].sum() == 0.0
    near = field[27:34, :].sum()
    assert near / field.sum() > 0.8
    assert stats["positive_output_cells"] > 0
