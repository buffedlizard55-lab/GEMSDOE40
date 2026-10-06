# H4-A negative result: Euler SI=0 depth-cluster contact network

**Date:** 2026-10-06 (UTC) · **Arm:** H4-A (rank 1 of the H4 slate) · **Status:** CLOSED — not promotable, no submission made.

## 1. What was built

The full H4-A pipeline (`scripts/run_depth_cluster_v2.py`, `src/gemsdoe40/depthcluster.py`,
`src/gemsdoe40/emit.py`) implements the arm exactly as pre-registered in
`docs/research/hypotheses.md`:

1. sliding-window Euler deconvolution at structural index N = 0 (fault-like contact) on
   challenge bands `rtp` (2) and `iso_grav_anom` (13), windows 8/12/16 px, with the gravity
   vertical derivative derived in the Fourier domain from the supplied field;
2. a 3-D depth-consistency gate (shallow × well-fitted × locally tight × mutually
   depth-consistent), magnetic and gravity solutions boosted where they agree in (x, y, z);
3. cross-field concordance pairing, weighted KDE splat to a continuous field `E`;
4. multi-scale Hessian lineament response `L` (trace geometry, not a point cloud);
5. surface-expression concurrence `C` from the 3DEP 1 m scarp-feature stack and band
   `det_elev_slope`;
6. H4-B mapping-gap weighting `G` from the official USGS QFFD scale split;
7. variants A0 = `E`, A1 = `L·(0.35 + 0.65C)`, A2 = A1·`G`;
8. 3 px non-maximum-suppressed dotted emission with the catalogue plus a 2 px buffer
   excluded, budget chosen from the measured precision ladder;
9. GeoTIFF twins (`nan` outside / `zeros` outside) plus zip, format re-read, and the
   registered uniqueness audit.

Solve counts: **120,099** magnetic and **49,037** gravity raw solutions → **83,930 / 28,672**
after the depth gate → **14,461** cross-field concordant pairs → KDE field with
**2,255,615** positive cells (43.6 % of the 5,167,373-cell footprint).

## 2. Result — the detector is beaten by chance on its own registered comparator

Truth and domain are the registered catalogue-block instrument: `data/labels.tif` (the
organizers' existing-fault raster, 60,988 positive cells), 4 rows × 6 columns of contiguous
blocks with a three-cell guard, pooled official distance-weighted Tversky index
(alpha = 0.2, beta = 0.8, R = 3 px).

| field | pooled DTI @ 90,000 dots | note |
| --- | --- | --- |
| A0 `euler_only` | 0.0172 | |
| A1 `euler_x_terrain` | 0.0169 | selected (tie rule: top two within 0.005) |
| A2 `euler_x_terrain_x_gap` | 0.0168 | |
| **same-mass random control** | **0.1257** | 90,000 uniform in-footprint cells |
| gradient-top-K baseline | 0.0365 | ranking by the incumbent gradient field |
| contiguous (no NMS) top-90k extraction | 0.0372 | isolates the emitter from the field |

The candidate scores **7.4× below the same-mass random control**. 2-fold held-out block
selection (choose the budget on one half, report on the other) gives 0.0106–0.0107, so the
result is not a budget-selection artefact. The pre-registered tie-break rule fired, so the
terrain-concurrence variant A1 was emitted; the emitted artifact is
`docs/downloads/gems40-h4a-euler-sicontact-depthcluster-20261006T012528Z-570e6300-nan.tif`.

## 3. Diagnosis — the field is nearly uniform where it must be selective

| measurement | value |
| --- | --- |
| AUC of the Euler field against the catalogue (0 / 1 / 2 px tolerance) | 0.534 / 0.538 / 0.541 |
| mean distance to nearest catalogue pixel, top-90k cells | 27.1 px |
| mean distance to nearest catalogue pixel, all positive cells | 28.1 px |
| mean distance to nearest catalogue pixel, random footprint cells | 30.8 px |
| KDE support | 2,255,615 px (43.6 % of the footprint) |

The 99.5-quantile normalisation leaves 44 % of the footprint positive, so the emitter is
choosing among essentially equivalent cells: the field separates catalogue-adjacent terrain
from the rest by less than 4 px of mean distance. This is a property of the *field*, not of
the emitter — the contiguous extraction without any NMS also lands at 0.0372. The failure
mode is therefore upstream: neither the depth gate nor the cross-field concordance nor the
lineament transform produced a spatially selective contact field, which is consistent with
the earlier layer-skill diagnostic in this repository (`docs/data/layer_skill_diagnostic.json`:
all 27 layers 0.42 ≤ AUC ≤ 0.60 against the catalogue).

## 4. What this does *not* show

The catalogue-block instrument measures skill against a *published, coarse-scale* fault
compilation (95.4 % of catalogue cells lie on USGS QFFD coarse-scale traces), while the
hidden target is, by the organizers' own definition, fault pixels **not** in USGS/INGENIOUS.
A candidate can therefore fail this comparator for reasons that have nothing to do with the
private metric. The correct reading of this result is: *H4-A cannot be shown to work on the
only instrument this repository registered before the arm was implemented*, so it is not
promotable — not that a similar physical approach must fail on the hidden set.

## 5. Attribution integrity (why no score-improving claim is made here)

An intermediate analysis in this session ranked the ten cached prior GeoTIFFs by
*group-reported* leaderboard scores and found strong correlations with emitted-dot
properties (distance to catalogue rho +0.988, detrended elevation +0.988, emitted pixel
count −0.976). That analysis is **withdrawn as evidence**. A leaderboard row establishes an
account-level public score only; it does not identify a raster unless the organizers link
the submission file, and `docs/data/feed-20261005.json` already prohibits attributing the
0.2778 row to a specific file. Both sides of those correlations are therefore unauthenticated.
The correlations are retained only in `docs/data/unauthenticated_attribution_audit.json` as a
documented record of the prohibited method, and no emission weight, threshold or design
decision in this repository may be derived from them. The draft design that had been built
on those weights was abandoned before any raster was written.

## 6. Disposition

- H4-A: **closed, not promotable.** Artifact moved out of `docs/downloads/` so it cannot be
  mistaken for a shippable candidate; retained under `work/withdrawn/` (gitignored).
- H4-B/G weighting is retained as a registered covariate only, not as a scoring claim.
- H4-C (radiometric alteration), H4-D (blocked by sandbox egress), H4-E (blind
  basin-interior) remain registered and unimplemented.
- The frozen promotion gate cannot be evaluated in this sandbox: its proxy raster
  (SHA-256 `26d142c4…`) lived in `/tmp/proxy-data/`, which does not persist between turns,
  and the surviving SGMC raster (`643cbe99…`) is the circular one used to build the frozen
  incumbent. Until the pinned proxy or an organizer-linked instrument is available, no arm
  can be promoted.

## 7. Audit status of the frozen artifact

The in-pipeline uniqueness pass did not complete inside the sandbox time budget, so the
registered checks were run separately on the frozen artifact
(`scripts/audit_h4a_artifact.py` → `docs/data/h4a_uniqueness_audit.json`):

- format re-read: `ok_shape`, `ok_range`, `ok_single_band`, `ok_float32` all true; EPSG:32611,
  3730 × 3292, float32, NaN outside the footprint, values in [0, 1], 90,000 positive px;
- uniqueness: `is_new = true` against all 278 same-grid cached priors (289 blobs minus
  grid mismatches), worst |Pearson| 0.0000 and worst positive-pixel Jaccard 0.0346, both far
  inside the registered thresholds (0.85 / 0.50).

Being unique and well-formed does not make it a good candidate: it is withdrawn on the
comparator result in section 2.
