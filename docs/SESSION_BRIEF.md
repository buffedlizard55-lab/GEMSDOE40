# Standing brief (read at the start of every session)

This is the operational extract of the user prompt. The full competition
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
