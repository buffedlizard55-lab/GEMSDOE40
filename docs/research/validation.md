# Append-only validation log — GEMSDOE40

This file records outcomes against [`hypotheses.md`](hypotheses.md). It does not rewrite a registered method, threshold, comparator, or decision after a result. Full evidence is linked by experiment; scores are local owner-derived SGMC proxy values, not organizer scores.

## 2026-10-05 UTC — H1 feasibility stop

- **Registered method:** Rank 1, joint magnetic–gravity Euler depth-consensus KDE; fixed 45,962-cell output budget.
- **Result:** 2,305 accepted magnetic solutions; 0 gravity solutions; 0 pairs. The budget could not be met.
- **Action:** abort; do not relax or backfill. No candidate TIFF and no holdout score.
- **Evidence:** [`feasibility-20261005.md`](../reports/feasibility-20261005.md), [`feasibility.json`](../data/feasibility.json).

## 2026-10-05 UTC — strict H2 feasibility stop

- **Registered method:** four-height TMI Euler persistence with fixed 45,962-cell budget.
- **Result:** 162 persistent anchors, 17,425 positive KDE cells; budget not met.
- **Action:** strict H2 aborted rather than backfill. No strict-H2 candidate or holdout score.
- **Evidence:** [`feasibility-20261005.md`](../reports/feasibility-20261005.md), [`feasibility.json`](../data/feasibility.json).

## 2026-10-05 UTC — H2-B proxy holdout failure

- **Registration:** after strict H2 failed its support gate and before any proxy score was read, a separate append-only H2-B rule was registered: emit every positive KDE-support cell and no other cell. It is not a rewrite or success of strict H2.
- **Candidate:** [`gemsdoe40-h2b-tmi-euler-natural-support-20261005.tif`](../downloads/gemsdoe40-h2b-tmi-euler-natural-support-20261005.tif), SHA-256 `02486eaa491da2d8ebc6b1cd53ed2bdc8e995f316be1110a64eec484779b5a68`.
- **Initial run:** Format PASS. Uniqueness PASS against the then-current 268 same-grid prior rasters. A later completeness addendum expanded the audit to 277 comparisons; the candidate still passes with 0 near-duplicates.
- **Proxy score:** 0.0063132852 pooled, 0/24 blocks won. Frozen all-prior incumbent 0.8359066541; delta −0.8295933689. The gate fails; action is **HOLD; DO NOT SUBMIT**.
- **No-go consequence:** no weekly slot was used, no upload was made, and there is no DrivenData score.
- **Proxy caveat:** the owner-derived SGMC mirror is not independently rebuilt truth. After score generation, the incumbent was identified as `gapfinder-v2-sgmc-gap`, a pinned prior whose originating project describes it as an SGMC-gap-only product. The all-prior comparator and all scores remain unchanged as preregistered. The comparator's score is circular on this proxy and is not an independent generalization estimate; this caveat was appended without altering the gate.
- **Evidence:** [`validation-h2b-20261005.md`](../reports/validation-h2b-20261005.md), [`validation.json`](../data/validation-h2b-20261005.json), [`format_receipt.json`](../data/format_receipt.json), [`uniqueness_audit.json`](../data/uniqueness_audit.json), [`source-register.md`](source-register.md).

## State after the initial H2-B run

H1 and strict H2 were feasibility failures; H2-B was a proxy promotion failure; H3 remained blocked pending a usable-plane and spatial-coverage audit. The later prior-corpus completeness addendum below supersedes the initial prior-count/uniqueness number but does not change the H2-B score or no-go decision. No candidate is eligible for a weekly submission slot. Future work must be separately registered and must not use these proxy outcomes to retune a supposedly preregistered result.

## 2026-10-05 UTC — H3 ComCat source-coverage audit (not a model score)

- The locked official ComCat count query returned 34,716 events advertising a focal-mechanism product; 238 had magnitude ≥4.0.
- One event detail (`nc10085763`, M6.4) exposed two finite nodal planes. This is an existence check, not a representative sample.
- The usable nodal-plane count across the full query and the number of registered 24 blocks containing usable mechanisms remain unverified. The ≥50 usable / ≥8 blocks gate is therefore **not passed**. H3 was not implemented.
- The official web research tool could retrieve count/detail pages; direct sandbox Python HTTPS was closed during TLS setup, so a reproducible local bulk-field audit could not be run.
- Evidence: [`h3_source_audit.json`](../data/h3_source_audit.json). No proxy, challenge labels, or leaderboard scores were used in this source audit.

## 2026-10-05 UTC — prior-TIFF corpus completeness addendum

- During branch/main reconciliation, 11 TIFF path observations in pinned `origin/main` commit [`0e4a795d467f16555443ee0446875d7c32f51aa5`](https://github.com/buffedlizard55-lab/GEMSDOE40/commit/0e4a795d467f16555443ee0446875d7c32f51aa5) were found outside the first inventory. They represent 9 new unique Git blobs. The inventory now has 279 unique blobs, with 277 exact-grid outputs scored and 167 passing the strict local format gate.
- H2-B was re-audited against all 277 exact-grid prior outputs. Uniqueness remains **PASS** (0 near-duplicates); the candidate's 0.0063132852 proxy DTI, 0/24 block wins, frozen incumbent (0.8359066541), and **HOLD; DO NOT SUBMIT** decision are unchanged.
- The mainline H40-4 multi-scale artifact has 97,654 positive cells. It scored 0.1247145885 on the owner-derived SGMC proxy and won 0/24 blocks versus the frozen incumbent; it fails the strict local format gate because it lacks a NaN nodata tag. The prior main project’s LM-calibrated 0.34835 and 4/4 folds are from a different instrument and are not comparable. Neither result is an organizer score.
- The unchanged incumbent is explicitly SGMC-gap-derived, so its result remains circular on this proxy and is not an independent generalization estimate. No hypothesis, threshold, or submission decision was retuned.
- Evidence: [current H2-B validation report](../reports/validation-h2b-20261005.md), [expanded validation JSON](../data/validation-h2b-20261005.json), [inventory](../data/prior_raster_inventory.json), [metadata audit](../data/prior_raster_metadata_audit.json), [expanded uniqueness audit](../data/uniqueness_audit.json), and [dated project/leaderboard feed](../data/feed-20261005.json). The original first-pass report and 268-comparison audit remain available alongside the current report.


## 2026-10-06 — appended H4 outcome (previous records unchanged)

The offset-aware, rank-adaptive contact Euler + cross-window depth-KDE run completed on CPU. It produced an independently new, exact-format TIFF and 46,656 QC-passing depth solutions (20,617 clustered). Its SHA-256 and its cloud's SHA-256 were independently reproduced byte for byte.

A final reachable-history inventory covers 55 public repositories, 333 branch heads, 1,835 commits and 336 byte-verified TIFF-named artifacts. 333 are full-size arrays; technical fixtures/header-only artifacts and two partial-NaN storage variants are explicitly adjudicated, not silently omitted or fabricated. Every complete counterpart is separately compared. No near-duplicate under the unchanged H4 thresholds.

**HOLD — DO NOT SUBMIT.** H4 proxy DTI is 0.00690153557133729; H33-B2 is 0.09155026618057825; the frozen SGMC-derived incumbent is 0.8359066540883113. H4 wins no truth-bearing block against either. The historical 18-win gate is also mathematically infeasible with only 16 truth-bearing blocks, and its incumbent is circular. We preserve, disclose and do not retroactively relax either defect. No organizer score or competition submission exists for H4.

[Full result](../reports/h4-results-20261006.md) · [Exact block scores](../data/h4-validation.json) · [Full raw-output audit](../data/h4-uniqueness.json) · [Three-pass review](../reports/review-20261006.md). Next: independently register a feasible whole-fault-system validation protocol before any new model search.
