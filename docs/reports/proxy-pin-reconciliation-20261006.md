# Proxy-pin reconciliation — 2026-10-06 UTC

## Decision

The two hashes are not competing spellings for the same file. They are byte-distinct, same-grid SGMC-derived owner mirrors from different repository paths. Preserve them as **separate validation-instrument versions**. Historical results tied to `26d142…` remain historical; all new candidate validation in this run uses `643cbe…`. Neither is organizer ground truth.

## Pinned source bytes

| Instrument | SHA-256 | bytes | pinned owner source | Git blob |
|---|---|---:|---|---|
| Current `643cbe` | `643cbe992ef4ba37588fb469163ed8291e3ceb23d6c1f78a3cfaa462430c2da0` | 198,602 | [GEMSDOE24 `data/external/derived_sgmc_faults_100m_u8.tif`](https://github.com/buffedlizard55-lab/GEMSDOE24/blob/07345ea0604953d7efb858d9cfbc21e20c7aca0b/data/external/derived_sgmc_faults_100m_u8.tif) | `c794f6ee8dd278f9d836576dcb5d02239da1c4f7` |
| Historical H2-B `26d142` | `26d142c4c93282cd94f6950ab96f22aeff59fbbea523d43d662e76fa1b161b5c` | 213,034 | [GEMSDOE30 `data/external/derived_sgmc_faults_100m_u8.tif`](https://github.com/buffedlizard55-lab/GEMSDOE30/blob/1f9ac110d5f1a4fac7957fe814267a0be964dd20/data/external/derived_sgmc_faults_100m_u8.tif) | `2e0aeaaee91791e7a2212d5f8133e2b26650f75f` |

Each downloaded byte stream was rechecked against the stated Git blob SHA-1 and SHA-256. The prior H2-B loader's `26d142…` expectation therefore has an identifiable upstream source; it is not silently changed to the current bytes.

## Raster comparison

Both are one-band uint8 rasters on the exact sample grid: 3,730 × 3,292, EPSG:32611, transform `(100, 0, 243350, 0, -100, 4508550)`. The old file has a `255` nodata tag but no 255-valued cells; the current file has no nodata tag. Neither metadata difference changes the in-footprint binary values at a cell.

- Historical positives: 82,151 cells across the full stored grid.
- Current positives: 83,593 cells across the full stored grid.
- Different cells: 1,450.
- Pearson correlation over the stored arrays: 0.991230.
- Positive-support Jaccard: 0.982655.
- Exact off-catalogue proxy positives after masking the exact known-label pixels: 78,433 historical; 79,615 current.

These are similar versions, but not identical labels. The 1,450 changed cells are enough to move DTI and block-level results; pin each result to its hash.

## Same-candidate sensitivity check

The public H2-B TIFF was fetched from the pinned current GEMSDOE40 `main` tree, then rescored without rebuilding it. The public H40-4 TIFF was checked the same way. This separates proxy-version sensitivity from candidate-construction changes.

| Existing output | Candidate SHA-256 | `26d142…` pooled proxy DTI | `643cbe…` pooled proxy DTI | Evidence boundary |
|---|---|---:|---:|---|
| H2-B TMI Euler natural support | `02486eaa491da2d8ebc6b1cd53ed2bdc8e995f316be1110a64eec484779b5a68` | 0.006313285 | 0.006223113 | Both are proxy scores, not organizer scores; H2-B remains failed. |
| H40-4 multiscale Euler augmented | `33d893193f6eb009165b5c62bcef1c7c945531ed19ddb2c621150a52f0c16835` | 0.124714589 | 0.125939967 | Both are proxy scores; the TIFF lacks NaN nodata and is not a format-eligible submission. |

The H2-B score `0.006313` is now byte-reproducible on its named historical instrument. It must not be restated as the current-proxy score. The H40-4 delta likewise cannot be compared across versions as though truth were unchanged.

## Reproducibility

- Current bytes: `data/external/derived_sgmc_faults_100m_u8.tif` (ignored local input; hash in `docs/data/input_manifest.json`).
- Archived H2-B bytes: `data/external/archive/derived_sgmc_faults_100m_u8_h2b_26d142.tif` (ignored local input; hash in `docs/data/input_manifest.json`).
- Historical manifest preserved as [`input_manifest-20261005-h2b-legacy.json`](../data/input_manifest-20261005-h2b-legacy.json).
- Candidate rasters rescored from public GitHub paths are in the refreshed prior-raster inventory and are not copied or edited.

## Remaining limit

The owner-derived proxy lineage is not independently reconstructed from the official USGS/SGMC product and is not equivalent to the hidden public or private competition labels. It is suitable only for an explicitly named local diagnostic. New H7 evaluation uses the current `643cbe…` proxy, with all proxy and circularity caveats retained.
