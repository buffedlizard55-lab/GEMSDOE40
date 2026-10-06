# Historical operational extract (pre-H2-B branch)

> **Superseded as current instructions.** This main-branch note predates the H2-B registration, expanded 279-blob prior audit, and current no-go decision. In particular, its generic advice about an unset nodata tag/zeros outside and its proposed hidden-set-matching instrument were not established as the current organizer contract or validated truth source. Use the top-level [`README`](../README.md), [`docs/prompt.md`](prompt.md), [`format_receipt.json`](data/format_receipt.json), and append-only [`hypotheses.md`](research/hypotheses.md) as the current record. H2-B remains **HOLD — DO NOT SUBMIT**.

This is the prior operational extract of the user prompt. The full competition
problem is at https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/

## Maximize P(Win). Own the Outcome.

- Every weekly submission slot is an experiment, not a lottery ticket.
- Do not copy a previous GEMSDOE TIF unless it is for learning. The shipped
  file must be unique (Pearson < 0.85 and Jaccard < 0.50 vs every prior).
- Do not spend a slot on an idea that has not beaten the current holdout
  best *on the instrument that matches the hidden set*. Catalogue-proxy DTI
  is the wrong objective for off-catalogue truth.
- Predicted values must be in [0, 1]. Primary download = all-finite,
  nodata unset, zeros outside the footprint. That is the fix for the portal
  error `Predicted values must be in range [0, 1]`.
- The download must be the first thing on the site.
- No hallucinations. Every number traces to a fetched source or a file
  hashed in this repository.

## This session's unique method

Euler deconvolution depth-clustering (Reid, Allsop, Granser, Millett &
Somerton, Geophysics 1990): SI = 0 (fault-like contact) on RTP magnetics
and isostatic gravity → cloud of depth-labelled solutions → weighted KDE
(tight, shallow, mutually-consistent clusters score higher) → [0, 1]
float32 GeoTIFF, EPSG:32611, 100 m, 3730 × 3292.

## Grid (measured, not assumed)

- sample_submission.tif sha256 `2176d08e485aa2cd2860ce8df539db4faf4d76163b38a4dd8c30a40454d35cbc`
- 3730 × 3292, EPSG:32611, transform (100, 0, 243350, 0, −100, 4508550)
- footprint 5,167,373 finite cells
