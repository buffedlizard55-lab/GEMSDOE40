# H7 candidate validation — 2026-10-06 UTC

## Decision

**HOLD — RESEARCH ONLY — DO NOT SUBMIT. No weekly slot was used.** The candidate is a valid and pixel-distinct TIFF, but its 24-block score against the current owner-derived SGMC proxy is negligible and loses every nonempty comparison block to the frozen highest-prior baseline. This is a local proxy diagnostic, not a DrivenData score and not evidence about private/final labels.

## Candidate and construction

- **Hypothesis:** H7, gravity-gradient-weighted magnetic Euler source cloud. It was registered before H7 implementation and built without reading either proxy raster.
- **Potential-field Euler solve:** 3,237 accepted `rtp` SI=0 source solutions under the fixed H4 magnetic gates; 430 positive-weight clustered 3-D voxels. Gravity `iso_grav_anom` contributes only the preregistered soft horizontal-gradient multiplier `[0.5, 1]`; no gravity Euler source is claimed.
- **Raster:** [`gemsdoe40-h7-rtp-euler-gravity-context-3d-kde-20261006-998f660f.tif`](../downloads/gemsdoe40-h7-rtp-euler-gravity-context-3d-kde-20261006-998f660f.tif); 28,080 positive in-footprint cells.
- **Raw SHA-256:** `e229aa9af26018bc80b0cca880d484a5b8362c26ad79ca362423953c32e66597`.
- **Canonical pixel SHA-256:** `998f660fcbf3b1a31600bd1e5e56a4d3b4c32d57c8fd8af36b3af793226a7d7c`.
- **Format:** PASS — one-band `float32`, 3,292 × 3,730, EPSG:32611, exact sample transform, 100 m north-up, 5,167,373 in-footprint cells, finite values `[0,1]`, and NaN in all 7,111,787 outside cells.

## Spatially blocked proxy result

The frozen 4-row × 6-column holdout with three-cell internal guards scored the candidate against the current mirror pinned to SHA-256 `643cbe992ef4ba37588fb469163ed8291e3ceb23d6c1f78a3cfaa462430c2da0`. Only exact `labels == 1` pixels were excluded. The current proxy is an owner-derived mirror, not authenticated organizer delivery or independently reconstructed ground truth.

| Comparator | Pooled proxy DTI | 24-block result | Interpretation |
|---|---:|---:|---|
| **H7 candidate** | **0.0002650893** | **0 wins** | Research diagnostic only |
| Highest prior, all outputs | 0.8257222257 | H7 delta −0.8254571364; 0 wins | `gapfinder-v2-sgmc-gap`, SHA-256 `7251c22b…e3097`; explicitly SGMC-gap-derived and circular on this proxy. Retained in the preregistered conservative gate. |
| Highest non-circular prior | 0.7618701028 | H7 loses | Geologic-map-faults-gap output, blob `4a0f0af3…`; local proxy score only, not organizer truth. |

The worst H7 block delta versus the all-prior incumbent was −0.9421647993 (registered limit −0.005). The score therefore fails decisively; no tuning or threshold change is authorized. The account-level public leaderboard and another owner's projected score are not used in this calculation.

## Prior-output completeness and uniqueness

The dated refresh scanned 55 public sibling-repository heads and 501 TIFF paths; every path was classified. The inventory has 306 unique Git blobs across 541 pinned locations, and all 306 were fetched and SHA-verified. The score scan compared H7 with **304 exact-template single-band prior rasters**; two other TIFFs were not on that grid/band. Of the 304, 183 also pass this project's strict local format filter.

H7 passed the preregistered uniqueness rule against all 304 comparable blobs: **0 near duplicates**, 0 raw-hash matches, and no canonical pixel match. The audit used a fixed 45,962-cell top-support budget and the frozen thresholds (absolute Pearson ≥0.995, nonzero-support Jaccard ≥0.90, or top-budget containment ≥0.90). The closest output by the audit ordering is the prior H2-B file (Pearson 0.052925, support Jaccard 0.052577, top-budget containment 0.130445). Two non-grid TIFFs were not pixel-compared but remain included in the complete, byte-verified inventory.

## Reproducible receipts

- Full block metrics, incumbent attribution, and gate decision: [`validation-h7-20261006.json`](../data/validation-h7-20261006.json).
- Scores for all 304 same-grid priors: [`prior-holdout-scores-h7-20261006.json`](../data/prior-holdout-scores-h7-20261006.json).
- Raw/canonical pixel uniqueness comparisons: [`uniqueness-audit-h7-20261006.json`](../data/uniqueness-audit-h7-20261006.json).
- Construction and independent format receipt: [`h7-build-20261006.json`](../data/h7-build-20261006.json) and [`format-receipt-h7-20261006.json`](../data/format-receipt-h7-20261006.json).
- Prior corpus inventory and refresh method: [`prior_raster_inventory-20261006.json`](../data/prior_raster_inventory-20261006.json) and [inventory refresh report](prior-raster-refresh-20261006.md).
- Preregistration, including the fixed top-budget and no-slot rule: [`preregistered-hypotheses-20261006.md`](../research/preregistered-hypotheses-20261006.md).

## Limitations and next steps

1. The proxy mirror's provenance is not independently established, and the highest all-prior output is circular on it. H7's poor result is still a strong reason not to spend a slot, but proxy performance is not private-label performance.
2. This is a geologic-contact signal, not a direct fault detector. Lithologic contacts, alteration, basin fill, and intrusions can create magnetic/gravity gradients without faults or geothermal productivity.
3. The public prior inventory cannot prove uniqueness against private/unpublished artifacts or identify which TIFF earned any account-level leaderboard score.
4. Continue only with new preregistered, physically motivated candidates and independently obtainable labels/holdouts. Re-running or tuning H7 against this proxy is not allowed. Any future candidate must beat the then-current incumbent and pass fresh format, cache-completeness, and uniqueness checks before human review. No file is cleared today.

**Organizer score:** none. **Upload:** none. **Weekly feedback slot:** unused.
