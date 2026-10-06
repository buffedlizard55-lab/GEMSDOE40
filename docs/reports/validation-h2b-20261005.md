# H2-B proxy holdout and artifact audit — expanded prior corpus, 2026-10-05 UTC

## Decision

> **HOLD — DO NOT SUBMIT.** H2-B passes the local GeoTIFF contract and remains unique against the expanded prior-output corpus, but decisively fails the frozen owner-derived SGMC-proxy promotion gate. No weekly submission slot was used, no upload was made, and there is no organizer score for this file.

## Candidate and exact artifact

- Hypothesis/output rule: append-only H2-B — four-height TMI Euler persistence KDE; emit all and only positive in-footprint KDE support outside exact known-label pixels, at binary 1.0. No top-K cutoff or backfill.
- File: [`gemsdoe40-h2b-tmi-euler-natural-support-20261005.tif`](../downloads/gemsdoe40-h2b-tmi-euler-natural-support-20261005.tif)
- SHA-256: `02486eaa491da2d8ebc6b1cd53ed2bdc8e995f316be1110a64eec484779b5a68`
- Canonical in-footprint pixels SHA-256: `a1af728d9fc0e3806b4e6a7d729e925f7c3941c5c076e4cc5d5be987317701d4`
- Positive cells: 17,425 / 5,167,373 (0.337%).
- Local format: **PASS** — one-band float32, 3,730 × 3,292, EPSG:32611, 100 m, exact sample transform/footprint, finite values in `[0,1]` inside, NaN outside. This is a repository-side check, not an organizer receipt. See the [format receipt](../data/format_receipt.json).
- Current uniqueness: **PASS** — 277 exact-sample-grid prior blobs compared; zero near-duplicates under the frozen raw/canonical hash, correlation, support-Jaccard, and top-budget containment rules. See the [expanded uniqueness audit](../data/uniqueness_audit.json).

## Prior-inventory completeness addendum

The initial H2-B run used the 270-blob inventory then present in this branch (268 exact-grid comparisons). During reconciliation with the advanced `origin/main`, 11 additional TIFF path observations were found in pinned commit [`0e4a795d467f16555443ee0446875d7c32f51aa5`](https://github.com/buffedlizard55-lab/GEMSDOE40/commit/0e4a795d467f16555443ee0446875d7c32f51aa5); these reduce to nine unique Git blobs. They are now included in the [pinned 279-blob inventory](../data/prior_raster_inventory.json) and [metadata audit](../data/prior_raster_metadata_audit.json). The expanded corpus contains 277 exact-grid prior outputs, of which 167 pass this repository's strict local format eligibility check. The initial validation and uniqueness snapshots remain preserved as [`validation-h2b-initial-20261005.json`](../data/validation-h2b-initial-20261005.json) and [`uniqueness-audit-h2b-pre-mainline-addendum-20261005.json`](../data/uniqueness-audit-h2b-pre-mainline-addendum-20261005.json).

The nine newly inventoried outputs were scored on the same frozen local proxy. All won 0/24 blocks against the frozen incumbent; none changed the incumbent or H2-B's result. These are diagnostics, not organizer scores:

| Mainline prior output (unique blob) | Positive cells | Local format gate | Proxy DTI diagnostic | Wins vs incumbent |
|---|---:|---|---:|---:|
| `GEMSDOE40-submission-nan.tif` (`76f9bf7a…`) | 45,784 | PASS | 0.01454197 | 0 / 24 |
| `GEMSDOE40-submission.tif` (`033c13e3…`) | 45,784 | FAIL: no NaN nodata / outside mask | 0.01454197 | 0 / 24 |
| `gems40-euler-augmented-incumbent…-zeros.tif` (`363f229b…`) | 60,710 | FAIL: no NaN nodata / outside mask | 0.09726689 | 0 / 24 |
| `gems40-euler-augmented-incumbent….tif` (`73f3074f…`) | 60,710 | FAIL: no NaN nodata | 0.09726689 | 0 / 24 |
| `gems40-euler-multiscale-augmented-incumbent…-zeros.tif` (`08a606d6…`) | 97,654 | FAIL: no NaN nodata / outside mask | 0.12471459 | 0 / 24 |
| `gems40-euler-multiscale-augmented-incumbent….tif` (`422ac946…`) | 97,654 | FAIL: no NaN nodata | 0.12471459 | 0 / 24 |
| `gems40-euler-si0-depthcluster-crossfamily…-zeros.tif` (`36465326…`) | 20,000 | FAIL: no NaN nodata / outside mask | 0.01107606 | 0 / 24 |
| `gems40-euler-si0-depthcluster-crossfamily….tif` (`2542597d…`) | 20,000 | FAIL: no NaN nodata | 0.01107606 | 0 / 24 |
| `gems40-euler-si0-depthcluster-crossfamily-continuous…tif` (`d4639f1b…`) | 3,380,827 | FAIL: no NaN nodata | 0.03729877 | 0 / 24 |

`H40-4` refers to the 97,654-cell multi-scale output. The pre-merge main project separately reported LM-calibrated DTI 0.34835 versus 0.26792 for its previous site-best and 4/4 favorable folds. That was a different local instrument and is not comparable to this SGMC-proxy DTI; neither value is an organizer score. H40-4's file also lacks the NaN nodata tag required by this repository's strict local format check, so its proxy score is diagnostic only.

## Frozen proxy comparison

The truth instrument is `/tmp/proxy-data/derived_sgmc_faults_100m_u8.tif`, SHA-256 `26d142c4c93282cd94f6950ab96f22aeff59fbbea523d43d662e76fa1b161b5c`. It is an **owner-derived SGMC mirror**, not an independently reconstructed official raster, organizer truth, or hidden challenge labels. Exact known-label pixels are masked only at their exact locations; internal block boundaries have a three-cell guard. Scores below are local Distance-Weighted Tversky values using the competition's 300 m triangular kernel and `alpha=0.2`, `beta=0.8`.

| Measure | H2-B candidate | Frozen all-prior incumbent |
|---|---:|---:|
| Pooled proxy DTI | **0.0063132852** | **0.8359066541** |
| Weighted TP / FP / FN | 411.306 / 16,342.772 / 76,836.694 | 62,072.707 / 225.000 / 15,175.293 |
| Positive cells | 17,425 | 61,664 |
| Block wins | **0 / 24** | — |
| Candidate minus incumbent | **−0.8295933689** | — |
| Worst block delta | **−0.9421647993** | — |

The promotion gate requires a strict pooled win, at least 18/24 block wins, and no block loss worse than −0.005. H2-B fails all three. The [complete machine-readable report](../data/validation-h2b-20261005.json) contains all 277 prior scores, candidate block components, frozen gate, and mainline addendum.

### Incumbent-provenance caveat

The frozen rule selects the highest pooled score among comparable, locally format-eligible exact-grid prior outputs. Its winner remains `GEMSDOE3/docs/downloads/gapfinder-v2-sgmc-gap-20260925T045014Z-7251c22bb4.tif` (pinned Git blob `9380203f3cb5edb9d096614f8dcf453e8f141b8f`; SHA-256 `7251c22bb489198ae4709611c29cef3ac799d8942c2f1779086d5b59ae0e3097`). Its name and originating project materials describe it as an SGMC-gap-only output, built from USGS SGMC faults more than 300 m from supplied labels. Because the local proxy is itself an owner-derived SGMC product, the incumbent's 0.8359 is a **circular proxy ceiling**, not an independent estimate of challenge generalization. It remains in the all-prior comparison as preregistered; the limitation is disclosed and the comparator is not represented as independent evidence.

Provenance references: [pinned prior TIFF](https://raw.githubusercontent.com/buffedlizard55-lab/GEMSDOE3/512371de5d30162813803cdccbd4931ceeb065c7/docs/downloads/gapfinder-v2-sgmc-gap-20260925T045014Z-7251c22bb4.tif), [originating GEMSDOE3 project](https://github.com/buffedlizard55-lab/GEMSDOE3), [USGS SGMC metadata](https://mrdata.usgs.gov/geology/state/USGS_SGMC_Metadata.xml). The official source metadata verifies the USGS compilation context, not the provenance or correctness of the owner-derived proxy raster used here.

## Experiment record and evidence boundaries

- H1 stopped at feasibility: 2,305 accepted magnetic solutions, zero accepted gravity solutions, no pairs, no candidate and no holdout score.
- Strict H2 stopped at its fixed 45,962-cell budget: it had 17,425 positive KDE cells and was not backfilled.
- H2-B's separate all-positive-support rule was appended before its proxy score was read. Its failure does not rewrite or rescue strict H2.
- H3 remains blocked: the official ComCat count query found 34,716 events advertising focal-mechanism products (238 at magnitude ≥4); usable nodal-plane count and spread across the registered 24 blocks remain unverified. H3 was not implemented.
- No hidden labels were accessed. SGMC was not used to construct, rank, or threshold H2-B; it was read only by the registered local proxy validator.
- The public leaderboard snapshot is separate: top four 0.3262, 0.3222, 0.3220, 0.3195 on 2026-10-05. The account-level 0.2778 row for `extradr19` is not tied to a file. See the [dated feed](../data/feed-20261005.json).
- No submission, weekly slot, organizer response, or official candidate score exists. **Do not submit H2-B.**

Related evidence: [feasibility report](feasibility-20261005.html), [format receipt](../data/format_receipt.json), [expanded uniqueness audit](../data/uniqueness_audit.json), [prior inventory](../data/prior_raster_inventory.json), and [source register](../research/source-register.md).
