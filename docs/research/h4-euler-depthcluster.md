# H4 — Euler deconvolution depth-cluster submission (final record)

**Method under test:** 3-D Euler deconvolution (Reid, Allsop, Granser, Millett & Somerton 1990,
DOI [10.1190/1.1442774](https://doi.org/10.1190/1.1442774)) with the structural index of a
fault-like contact (SI = 0) on the reduced-to-pole magnetic grid and the isostatic-residual gravity
grid of `training_features.tif`, deconvolved at three upward-continuation heights and three window
sizes, kept only where solutions persist across heights, weighted by shallowness, mutual depth
consistency, cluster tightness and analytic-signal strength, and splatted into a kernel-density
estimate of **solution density** — not a thresholded gradient, curvature or edge magnitude.

**Shipped file (portal file):**
`docs/downloads/gemsdoe40-euler-line-ring-pruned-60000px-20261006T032708Z-8dafb186-zeros.tif`
· 274,414 B · SHA-256 `8dafb1860367e7d5b8530f924b2ff20efbdff32cfffda82647711f0dbb2e646d`
· 60,000 positive pixels, values exactly 0.0 / 1.0, `nodata=None`, float32, EPSG:32611,
3730 × 3292, transform `|100, 0, 243350| / |0, −100, 4508550|`.

**Receipt:** `docs/data/h4_submission_20261006T032708Z.json`.
**Byte audit:** `docs/data/h4_shipped_audit_zeros.json`, `…_nan.json` (both PASS).
**Gate:** `docs/data/h4_blocked_validation_shipped.json` (shipped file),
`docs/data/h4_blocked_validation.json` (rejected alternatives).

Nothing on this page is an organizer score.

---

## 0. Gate result — the file is NOT promoted

`scripts/validate_h4_blocked.py` scores candidate rasters on the repository's frozen **LM
instrument**: four spatially blocked quadrants, off-catalogue SGMC truth, prevalence-calibrated,
domain eroded 12 px, catalogue masked pixel-exactly. Stage 1 verifies the instrument itself against
the published (owner-reported) score ladder: it reproduces **5 of 6** published orderings when
calibrated (raw LM only 2 of 6), so it is informative for this family of fields.

| field | LM | LM-calibrated | folds better than the best prior | projected leaderboard DTI |
|---|---|---|---|---|
| reference `h33-2-b2` (owner-reported 0.2778) | 0.0726 | **0.3663** | — | 0.2778 (assumed pairing) |
| **shipped H4 (ring-pruned field, 60,000 px)** | 0.0790 | 0.2860 | **1 / 4** | 0.197 ± 0.021 |
| H4 reserve: flank band 2–6 px, 30,000 px | 0.0257 | 0.1168 | 0 / 4 | 0.028 ± 0.021 |
| H4 reserve: flank band + corroboration, 25,000 px | 0.0207 | 0.0968 | 0 / 4 | 0.008 ± 0.021 |

The frozen promotion rule is "beat the incumbent in all four blocked folds". The shipped file wins
**one** fold (NW 0.2824 vs 0.2446) and loses three, so **H4 is not promoted and the weekly
submission slot is not spent on this evidence**. The file is published because the project brief
mandates a unique Euler depth-clustering submission artifact and because the byte-level and
uniqueness audits are worth having on disk; the site states the gate result next to every download
link.

**Budget sensitivity on the same gate** (ring-pruned field, 2√2 px lattice):

| emission | 25,000 px | 40,000 px | 60,000 px (shipped) | 90,000 px | 140,000 px | 200,000 px |
|---|---|---|---|---|---|---|
| LM-calibrated | 0.1762 | 0.2337 | 0.2860 | 0.3323 | 0.3692 | 0.4013 |

At 200,000 px the LM-calibrated *mean* passes the reference (0.4013 > 0.3663 on the 2√2 lattice;
0.5585 on a 4 px lattice) but the fold-wise rule still fails (2 of 4 folds) and, more importantly,
LM's denominator adds `+0.8·tp_g` where the real metric subtracts `0.2·tp_g`, so LM systematically
under-penalises mass. Its monotone mass preference is therefore **not** used as evidence: the
60,000 px budget comes from the emission law fitted to the published score ladder
(`ρ(score, emitted px) = −0.907`, optimum ≈ 60,069 px), which is the only mass rule in this
project derived from real scores rather than from a surrogate.

## 1. What was built

| stage | setting |
|---|---|
| fields | `training_features.tif` band `rtp` (reduced-to-pole magnetics) and band `iso_grav_anom` (isostatic residual gravity), both sha256-pinned |
| continuation heights | 0, 500, 1,500 m upward continuation |
| windows | 6, 8, 12 px; stride 4 px; analytic-signal percentile 72; max relative standard error 0.22 |
| solutions | RTP 388,173 → 308,944 after cross-height persistence (median depth 437.6 m); gravity 68,102 → 36,783 (median depth 1,131.3 m) |
| per-solution weight | shallowness × cluster tightness (2 px neighbourhood count) × mutual depth consistency × analytic-signal strength |
| density | Gaussian KDE of solution density, σ = 1.6 px (2,566,449 positive pixels) |
| texture | line response (half-length 3 px) — the depth-cluster field is a lineament field, so the matched filter runs **along** it |
| catalogue ring | every pixel within **2 px** of a published-catalogue pixel is zeroed before emission (`--ring-prune-px 2.0`) |
| emission | greedy minimum-separation lattice, 2√2 px, **60,000 dots**, binary 1.0, zero elsewhere |

## 2. Why this emission

1. **The method is the brief's method.** Depth-clustering of Euler solutions, not a gradient
   threshold: the ranking surface is a density of depth-labelled solutions weighted so that tight,
   shallow, mutually-consistent clusters outrank scattered or deep ones.
2. **The ring-prune is evidence-anchored, not aesthetic.** In the 19-field prior-geometry study
   (`docs/data/prior_geometry.json`, n = 19): the programme's highest-scoring file kept **zero** mass
   within 2 px of the catalogue; every prior with mass inside 2 px scores below every prior that
   avoided it at matched mass; Spearman(owner score, fraction inside 2 px) = **−0.589**;
   Spearman(owner score, emitted pixels) = **−0.872**. The shipped file keeps 0 % inside 2 px,
   median distance 21.2 px, 14.9 % inside 6 px, 85.1 % beyond 6 px.
3. **It is the most distinct field this programme has produced.** Against all 19 prior rasters:
   maximum |Pearson| **0.0128** (`h18-3a-x-complexity`), maximum support Jaccard **0.0093**
   (`h18-3a-x-complexity`), 0 near-duplicates. Verdict NEW.
4. **It is the best of the seven H4 ranking surfaces on the frozen gate.** Ring-pruned 0.2860 > the
   raw KDE and the flank-gated variants (all ≤ 0.117 at their best budgets).

## 3. Independent byte audit

`scripts/audit_submission.py` re-reads the shipped bytes and re-checks every claim: identical shape,
transform, CRS and footprint to `sample_submission.tif`; single-band float32; every finite value in
`[0, 1]`; **no `nodata` sentinel** (the cause of the portal error "Predicted values must be in range
[0, 1]"); no positive value on a published-catalogue pixel; SHA-256 agreement with the receipt;
uniqueness against all 19 priors. Both twins PASS; exit status is non-zero on any failure, so this is
CI-able. Results: `docs/data/h4_shipped_audit_zeros.json`, `docs/data/h4_shipped_audit_nan.json`.

**Twin convention.** `-zeros` is the portal file (finite everywhere, values 0/1 only).
`-nan` is the sample-convention twin (NaN outside the valid footprint) for reviewers who prefer the
sample's convention. They carry the same positive support.

## 4. What is NOT claimed

* **No organizer score.** The public leaderboard exposes no filenames or hashes; every
  (file, score) pairing used anywhere in this repository — including the 0.2778 reference — is an
  **assumption**, not a receipt. The projected leaderboard DTI (0.197 ± 0.021) is an affine
  projection of a proxy instrument fitted to nine such assumptions (R² 0.913, RMSE 0.0107); it is a
  projection, not a measurement.
* **The modelled 0.1579 is not a score either.** It is the output of a two-parameter regression on
  (surrogate credit, surrogate precision) across 18 priors (RMSE 0.031, R² 0.807, LOO RMSE 0.0381);
  for this geometry it is an extrapolation, and §5 shows it disagrees with the frozen gate.
* **No claim that this beats anything.** It does not pass the four-fold rule; it is a unique,
  format-valid, audited research artifact.
* **The SGMC mirror is a proxy.** `data/external/sgmc/derived_sgmc_faults_100m_u8.tif`
  (sha256 `d569d553…700f`) is an owner-derived mirror of an SGMC-derived product, not independent
  challenge truth. An earlier repository record quoted a different hash (`26d142c4…b5c`) for the same
  filename; every number on this page is tied to `d569d553…700f`, and the conflict is unresolved
  because no USGS-side original is reachable from the build sandbox.
* **The metric's mask convention is read from a staff clarification, not from the rules PDF**:
  the known-fault mask is pixel-exact (no dilation). The shipped file does not depend on it (0 px
  inside 2 px of the catalogue), but the reconstructed scoring model used in §5 does.

## 4b. Mirror-version robustness (re-checked 2026-10-06)

The gate's truth layer is not stable: the same pinned repository path has yielded three different
hashes — `d569d553…700f` (receipts below), `26d142c4…b5c` (earlier record), `643cbe99…`
(re-fetched 2026-10-06). Re-running the gate against the newest mirror moves every absolute number
but not the outcome:

| mirror | reference `h33-2-b2` | shipped H4-line | folds | outcome |
|---|---|---|---|---|
| `d569d553…700f` | 0.3663 | 0.2860 | 1 / 4 | not promoted |
| `643cbe99…` | 0.2679 | 0.2085 | 1 / 4 | not promoted |

Absolute scores and the affine projection are mirror-specific and are never used to select a file;
only the invariant binary outcome is. Raw re-check:
`docs/data/h4_blocked_validation_recheck_20261006.json`.

## 5. Rejected alternatives, and the instrument disagreement that rejected them

The emitter's own surrogate model **preferred** the flank-band family: `euler_line_flank_gated`
(30,000 px, 4 px lattice) modelled 0.3410, `euler_line_flank_gated_corrob` (25,000 px) 0.3370 —
both above the 0.2824 modelled incumbent. The frozen gate rejects both (LM-cal 0.1168 and 0.0968;
0 of 4 folds, projected 0.028 / 0.008). Diagnosis, in order of evidential weight:

1. **No corpus analogue.** Every one of the 19 priors in the calibration spreads its mass across the
   footprint. The flank-band fields place **93–96 %** of their dots inside 2–6 px of the catalogue.
   No calibration point has that geometry, so the two-parameter model is extrapolating when it
   praises them; the gate, which computes a metric-shaped ratio for each field, is not.
2. **The SGMC surrogate rewards catalogue-hugging, the real ladder punishes it.** Only 12.8 % of the
   off-catalogue SGMC truth lies in the 2–6 px band, but that band is 3.9 % of the footprint, so the
   band is 3.3× enriched — a co-location artefact of two different agencies digitising the same
   structures. Across the corpus, the direct SGMC metric is **anti-correlated** with the published
   scores (Spearman **−0.921**). Any field that maximises SGMC credit is therefore buying a
   surrogate artefact.
3. **Both flanks of the family are therefore excluded from auto-selection.** `run_h4.py` now refuses
   to auto-select a surface with more than 50 % of its mass within 6 px of the catalogue unless a
   human forces it with `--surface`; the receipt records `emitted_distance_band_fractions`,
   the shipped geometry, and the gate reference for every emission.

The reserve files were deleted after measurement; their receipts
(`docs/data/h4_submission_20261006T032134Z.json`, `…T032157Z.json`) and the gate results
(`docs/data/h4_blocked_validation.json`) are retained as the audit trail.

## 6. Reproduce

```bash
PYTHONPATH=src python3 scripts/build_h4_clouds.py                 # Euler clouds (≈2 min)
PYTHONPATH=src python3 scripts/run_h4.py --stage field            # ranking surfaces (≈40 s)
PYTHONPATH=src python3 scripts/run_h4.py --stage emit \
    --surface euler_line_ring_pruned --spacing 2.8284 --budget 60000
PYTHONPATH=src python3 scripts/validate_h4_blocked.py --candidates <shipped-zeros.tif> \
    --out docs/data/h4_blocked_validation_shipped.json
PYTHONPATH=src python3 scripts/audit_submission.py <shipped-zeros.tif> \
    --receipt docs/data/h4_submission_<stamp>.json --out docs/data/h4_shipped_audit_zeros.json
```

The emission is deterministic: the same command and the same surfaces reproduce the same
`8dafb186…` digest.

Gate validation needs `data/existing_faults.tif`. It is absent from the public checkout; the
instrument reads `data/labels.tif` through a symlink shim, and the SGMC raster is read from
`data/external/sgmc/`. Both paths are configurable through `scripts/validate_h4_blocked.py --data`.

## 7. Limitations and remaining work

See [`docs/limitations.html`](../limitations.html) for the standing list. The ones that bind this
file: the missing fresh-DEM path (no network egress to USGS 3DEP from the build sandbox), the
absence of an independent non-mirrored truth layer, the unresolved SGMC hash conflict, and the fact
that no candidate in the programme has yet passed the four-fold rule.
