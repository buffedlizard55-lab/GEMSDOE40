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
- **Format:** PASS. **Uniqueness:** PASS against 268 same-grid prior rasters.
- **Proxy score:** 0.0063132852 pooled, 0/24 blocks won. Frozen all-prior incumbent 0.8359066541; delta −0.8295933689. The gate fails; action is **HOLD; DO NOT SUBMIT**.
- **No-go consequence:** no weekly slot was used, no upload was made, and there is no DrivenData score.
- **Proxy caveat:** the owner-derived SGMC mirror is not independently rebuilt truth. After score generation, the incumbent was identified as `gapfinder-v2-sgmc-gap`, a pinned prior whose originating project describes it as an SGMC-gap-only product. The all-prior comparator and all scores remain unchanged as preregistered. The comparator's score is circular on this proxy and is not an independent generalization estimate; this caveat was appended without altering the gate.
- **Evidence:** [`validation-h2b-20261005.md`](../reports/validation-h2b-20261005.md), [`validation.json`](../data/validation.json), [`format_receipt.json`](../data/format_receipt.json), [`uniqueness_audit.json`](../data/uniqueness_audit.json), [`source-register.md`](source-register.md).

## State after this log

H1 and strict H2 are feasibility failures; H2-B is a proxy promotion failure; H3 remains blocked pending a usable-plane and spatial-coverage audit. No candidate is eligible for a weekly submission slot. Future work must be separately registered and must not use these proxy outcomes to retune a supposedly preregistered result.

## 2026-10-05 UTC — H3 ComCat source-coverage audit (not a model score)

- The locked official ComCat count query returned 34,716 events advertising a focal-mechanism product; 238 had magnitude ≥4.0.
- One event detail (`nc10085763`, M6.4) exposed two finite nodal planes. This is an existence check, not a representative sample.
- The usable nodal-plane count across the full query and the number of registered 24 blocks containing usable mechanisms remain unverified. The ≥50 usable / ≥8 blocks gate is therefore **not passed**. H3 was not implemented.
- The official web research tool could retrieve count/detail pages; direct sandbox Python HTTPS was closed during TLS setup, so a reproducible local bulk-field audit could not be run.
- Evidence: [`h3_source_audit.json`](../data/h3_source_audit.json). No proxy, challenge labels, or leaderboard scores were used in this source audit.
