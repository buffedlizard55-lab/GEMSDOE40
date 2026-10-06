# Preregistered feasibility stops — 2026-10-05 UTC

## H1: joint magnetic–gravity Euler depth-consensus KDE

**Outcome: stopped at the frozen feasibility/emission gate. No candidate TIFF or holdout score.**

Using the fixed SI=0/10×10/stride-5 solve and preregistered residual, condition, horizontal-offset, depth, and pairing rules, H1 produced 2,305 accepted magnetic solutions and zero accepted gravity solutions. Therefore it had zero cross-field pairs, zero paired KDE support, and could not emit the fixed 45,962 positive cells. The runner aborted rather than relaxing criteria or backfilling.

Checks and inputs:

- Named `tmi_vg` validation passed its lock: Pearson `r = 0.9959886` and least-squares scale ratio `0.988864` against the reflect-padded Fourier downward-positive derivative.
- The common real-input mask had 5,164,312 cells; 60,988 exact known-label pixels were suppressed. No 300 m mask was used.
- The accepted magnetic cloud had p05/median/p95 depth about 238.59/551.84/1,161.82 m and median/p95 residual about 0.1645/0.1971.
- Gravity had no accepted solutions under the lock. This is a feasibility result, not a holdout result and not evidence that the method beats or loses to any baseline.

A diagnostic-only pass separately relaxed residual/depth/horizontal acceptance checks while retaining the condition-number gate. It found 202,453 solve windows for each field. Gravity depth median/p95 was about 15.14/95.61 m; residual median/p95 was about 0.7508/0.9604. This is recorded only to understand the failure; it did not generate a candidate, choose thresholds, or produce a score.

## H2: four-height TMI Euler persistence with fixed 45,962-cell budget

**Outcome: stopped at its fixed-budget gate. No strict-H2 candidate TIFF or proxy holdout score.**

The locked four-height pipeline produced 2,214 accepted height-0 solutions, 7,698 at 500 m, 4,512 at 1,000 m, and 1,629 at 2,000 m. After the locked 3-of-4 persistence match, 162 anchors remained (148 matches at three heights; 14 at four). The KDE had only 17,425 eligible positive cells, below the preregistered 45,962-cell emission budget. It aborted rather than backfilling.

Before reading a proxy score, a distinct append-only H2-B rule was registered to retain all and only the natural positive KDE support. H2-B is not a retroactive pass for strict H2. Its candidate was format- and uniqueness-audited, then scored; it failed and is marked HOLD in the [validation report](validation-h2b-20261005.md).

## Decision

Neither H1 nor strict H2 passed its own feasibility gate. H2-B did not pass its inherited promotion gate. **No weekly submission slot was used.** Full machine-readable counts and pinned input hashes: [`feasibility.json`](../data/feasibility.json). The H2-B candidate and all reports are research artifacts only; do not submit them.
