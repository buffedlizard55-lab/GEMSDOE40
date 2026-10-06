# Initial H2-B proxy holdout and artifact audit — 2026-10-05 UTC

> **Historical first-pass snapshot.** This report predates the current-main TIFF inventory addendum. It records 268 same-grid comparisons, not the final 277. The candidate result and decision did not change; use the [current expanded validation report](validation-h2b-20261005.md) and current [validation JSON](../data/validation-h2b-20261005.json) for the final audit.

## Decision

> **HOLD — DO NOT SUBMIT.** No weekly submission slot was used. The candidate passes the local GeoTIFF format checks and the frozen raw/canonical uniqueness checks, but it decisively fails the registered proxy promotion gate. It has no DrivenData score.

## Candidate

- Hypothesis/output rule: append-only H2-B — four-height TMI Euler persistence KDE; emit all and only positive in-footprint KDE support outside exact known-label pixels, at binary 1.0. No top-K cutoff or backfill.
- File: [`gemsdoe40-h2b-tmi-euler-natural-support-20261005.tif`](../downloads/gemsdoe40-h2b-tmi-euler-natural-support-20261005.tif)
- SHA-256: `02486eaa491da2d8ebc6b1cd53ed2bdc8e995f316be1110a64eec484779b5a68`
- Canonical in-footprint pixels SHA-256: `a1af728d9fc0e3806b4e6a7d729e925f7c3941c5c076e4cc5d5be987317701d4`
- Positive cells: 17,425 / 5,167,373 (0.337%).
- Format: **PASS** — one-band float32, 3,730 × 3,292, EPSG:32611, 100 m, exact sample transform and footprint, finite `[0,1]` inside, NaN outside. Full [format receipt](../data/format_receipt.json).
- Uniqueness: **PASS** — compared against 268 exact-sample-grid prior TIFF blobs. No raw/canonical hash match, `|r| ≥ 0.995`, support Jaccard ≥0.90, or top-support containment ≥0.90. Full [uniqueness audit](../data/uniqueness-audit-h2b-pre-mainline-addendum-20261005.json).

## Frozen proxy comparison

The truth instrument is `/tmp/proxy-data/derived_sgmc_faults_100m_u8.tif`, SHA-256 `26d142c4c93282cd94f6950ab96f22aeff59fbbea523d43d662e76fa1b161b5c`. It is an **owner-derived SGMC mirror**, not an independently reconstructed official raster, organizer truth, or hidden challenge labels. Exact known-label pixels are masked only at their exact locations; the internal block boundaries have a three-cell guard. Scores below are local Distance-Weighted Tversky values using the competition's 300 m triangular kernel and `alpha=0.2`, `beta=0.8`.

| Measure | H2-B candidate | Frozen all-prior incumbent |
|---|---:|---:|
| Pooled proxy DTI | **0.0063132852** | **0.8359066541** |
| Weighted TP / FP / FN | 411.306 / 16,342.772 / 76,836.694 | Not repeated here; see full prior score table in `validation.json` |
| Positive cells | 17,425 | 61,664 |
| Block wins | **0 / 24** | — |
| Candidate minus incumbent | **−0.8295933689** | — |
| Worst block delta | **−0.9421647993** | — |

The promotion gate requires a strict pooled win, at least 18/24 block wins, and no block loss worse than −0.005. H2-B fails all three. The per-block receipt is in [`validation.json`](../data/validation-h2b-initial-20261005.json).

### Incumbent-provenance caveat

The literal frozen rule selects the best pooled score from all comparable, format-eligible same-grid prior outputs. Its winner is `GEMSDOE3/docs/downloads/gapfinder-v2-sgmc-gap-20260925T045014Z-7251c22bb4.tif` (pinned Git blob `9380203f3cb5edb9d096614f8dcf453e8f141b8f`; SHA-256 `7251c22bb489198ae4709611c29cef3ac799d8942c2f1779086d5b59ae0e3097`). The artifact's name and originating project materials describe it as an SGMC-gap-only output, built from USGS SGMC faults more than 300 m from the supplied labels. Because the proxy truth is itself an owner-derived SGMC product, the incumbent's 0.8359 is a **circular proxy ceiling**, not an independent estimate of challenge generalization. This repository keeps it in the all-prior comparison exactly as preregistered, discloses the circularity, and does not relabel it as a fair independent benchmark. The candidate fails by a very large margin even against this conservative frozen comparison; no promotion follows.

Provenance references: [pinned prior TIFF](https://raw.githubusercontent.com/buffedlizard55-lab/GEMSDOE3/512371de5d30162813803cdccbd4931ceeb065c7/docs/downloads/gapfinder-v2-sgmc-gap-20260925T045014Z-7251c22bb4.tif), [originating GEMSDOE3 project](https://github.com/buffedlizard55-lab/GEMSDOE3), [USGS SGMC metadata](https://mrdata.usgs.gov/geology/state/USGS_SGMC_Metadata.xml). The official source metadata verifies the USGS compilation, not the provenance or correctness of the owner-derived proxy raster used in this run.

## What the experiment found

- The four accepted Euler clouds contained 2,214 height-0, 7,698 height-500 m, 4,512 height-1,000 m, and 1,629 height-2,000 m solutions.
- Only 162 height-0 anchors matched at three or four continuation heights (148 at three; 14 at four). Their weighted KDE supported 17,425 positive cells.
- The fixed H2 top-45,962 preregistered gate therefore stopped, rather than backfilling. The separate H2-B all-positive-support rule was appended to the register **before** H2-B's proxy score was read; it does not rewrite the failed H2 rule.
- Sparse support did not translate into proxy coverage: weighted TP was 411.3 against 77,248 scored holdout truth pixels, while the map paid 16,342.8 weighted FP. This is evidence against submitting this candidate under the frozen proxy criterion, not proof that the method has no geological value.

H1 separately stopped at feasibility with 2,305 accepted magnetic solutions, zero accepted gravity solutions, and no pair cloud; it was never holdout-scored. See [feasibility report](feasibility-20261005.md) and machine-readable [feasibility data](../data/feasibility.json).

## Evidence boundaries

- This is not an official score and does not estimate the private leaderboard without substantial proxy uncertainty.
- No hidden labels were accessed. SGMC was not used to construct, rank, or threshold the H2-B candidate; it was read only by the registered local proxy validator.
- The public leaderboard snapshot is separate: top four 0.3262, 0.3222, 0.3220, 0.3195 on 2026-10-05. The account-level 0.2778 row for `extradr19` is not tied to a file. See [`feed.json`](../data/feed.json).
- No upload, submission slot, or organizer response occurred. **Do not submit H2-B.**
