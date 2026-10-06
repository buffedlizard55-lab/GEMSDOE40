# Concurrent mainline integration — 6 October 2026

While contact-offset H4 was being implemented and independently reproduced, other sessions merged PRs #7 (H4-A), #8 (H7) and #9 (H40). This branch fetched and integrated main at `9987e40` without switching branches or deleting those contributions. Their source, tests, receipts, and committed TIFFs remain present.

## Distinct experiments, not aliases

| Name | Construction / record | Interpretation here |
|---|---|---|
| **Contact-offset H4** | TMI + dG/dz, SI=0 with free A, identifiable-rank solve, cross-window depth KDE | Current new research download. Native continuous field, no dotted emitter, no score-attribution calibration. HOLD. |
| **H4-A contact network** | Separate RTP/gravity/terrain/QFFD arm | Closed negative result on its own catalogue comparator. Its raw TIFF was moved to another session's ignored `work/withdrawn/`; no accessible raw file is invented or claimed as compared. |
| **H7 gravity-context** | Magnetic RTP Euler with a soft gravity-gradient multiplier | Retained with original bytes/receipts. Its registered proxy version differs from contact-offset H4. Not a dual-field Euler result and not an official score. HOLD. |
| **H40 depth-cluster / emission** | Separate cloud method and dotted-emission policy | Both committed TIFF twins are retained and included in the refreshed raw-output corpus. Its promotion estimates were withdrawn by that session. HOLD. |

“H7” in the contact-offset slate means **TMI/RTP representation stability**, still unimplemented in that slate. It is not the independently implemented H7 gravity-context arm. The registrations remain unchanged; this note resolves the naming collision rather than rewriting either record.

## Proxy versions must not be mixed

- Contact-offset H4 deliberately reproduces the frozen **26d142c4…** proxy, with incumbent DTI **0.8359066540883113**. That exact file was recovered from public GEMSDOE30 commit `1f9ac110…`, resolving the former missing-input blocker. It is still circular and the 18-win rule is infeasible with 16 truth-bearing blocks.
- H7/H40's retained manifest and `gems40.pins` describe the different **643cbe99…** profile. Their native scores must not be pooled with H4's as if the truth rasters were identical.
- The shared holdout helper now supports explicit version pins while preserving the historical default. H4's acquisition script independently verifies its frozen four inputs. A caller using another experiment must pass that experiment's proxy explicitly; silently replacing `data/external/derived_sgmc_faults_100m_u8.tif` with a different same-named file is forbidden.
- The inherited instrument is **retired for promotion**, not erased from history. No alternative local gate is adopted by this integration.

## Score-attribution irregularity

The incoming 16-pair instrument audit describes some scores as organizer-scored, while the parallel H4-A record explicitly retracts calibration using unauthenticated file-score associations. Both records are preserved. This session has no organizer file-level receipt authenticating those pairings. Their correlations and inferred hidden-label counts are therefore **conditional owner-attributed analyses, not independently verified official facts**, and are not used to choose H4's method, weights, threshold, or emission.

The proposed “achieve leave-one-out Spearman ≥0.8 on those pairs” rule is not accepted as a replacement for independent fault-system validation. Selecting a model to rank an unauthenticated, repeatedly reused score list can itself overfit.

## Integration safeguards

- Primary Pages content is generated from the H4 receipts. Secondary links preserve H7/H40 downloads with explicit no-upload warnings and the H4-A withdrawal record.
- The H7-only corpus refresher is preserved as `scripts/refresh_prior_inventory_h7.py`; the main refresher additionally inspects ZIP members and works with the reachable-history census.
- The entirely new upstream H40 runner is retained, not replaced with its older implementation. It requires deliberate research opt-in because it is not the current H4 release or an approved upload path.
- Upstream legacy uniqueness budget support is combined with exact-grid and empty/missing-cache fail-closed checks. Incoming tests are retained.
- The prior corpus is refreshed **before this branch's first push**, so its public-head snapshot includes newly merged H7/H40 and other newly accessible outputs, but does not accidentally compare H4 to its own just-published release.
- All model settings, H4 TIFF bytes and cloud bytes remain unchanged through integration. Source-only holdout-helper changes add explicit proxy-profile metadata, not a new scoring formula.

See the [current full audit](../data/h4-uniqueness.json), [three-pass ledger](review-20261006.md), [incoming README snapshot](../archive/README-main-9987e40-20261006.md), and [current evidence page](../evidence.html).
