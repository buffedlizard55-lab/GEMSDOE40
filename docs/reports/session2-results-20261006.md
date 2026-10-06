# Session 2 results — H8 trace-locked and H13 crest-binary Euler emissions — 6 October 2026

**Decision for both candidates: HOLD — DO NOT SUBMIT.** No weekly slot used. No organizer score.

## What was tested

The brief's required method — Euler deconvolution (Reid et al. 1990, SI = 0 fault-like contact)
producing a depth-labeled solution cloud, converted to a raster by kernel density of shallow,
mutually consistent solutions — was implemented in session 1 as **H4** and held (proxy DTI
0.0069). Session 1's own diagnosis of that failure: the KDE emission is spatially diffuse and its
arrangement does not resemble any ≥0.24 live artifact. Session 2 registers five new hypotheses
([H8–H12 slate](../research/h8-preregistration-20261006.md)) and implements the top two:

| Candidate | Construction | Frozen before scoring |
|---|---|---|
| **H8** | Frozen H4 cloud → per-family KDE → consensus field → structure-tensor orientation → orientation-NMS crests → along-strike depth consensus (±800 m × ±300 m corridor, ≥3 solutions, weighted depth SD ≤ 445 m) → continuous `C·γ` field, σ = 1 px, exact known pixels zeroed | yes |
| **H13** | Identical physics; emission restricted to binary 1.0 on the depth-consensus crest cells (metric-derived: linear value charging makes concentrated mass dominate diluted flanks) | yes (registered after H8's gate result, before any H13 scoring; addendum appended to the same preregistration) |

Both read only: the hash-pinned H4 cloud CSV (`6bed30b2…`, 46,656 records), the sample template,
and exact known labels. **No proxy truth, prior prediction, or leaderboard value enters
construction.**

## Gate results (frozen gates G1–G5)

| Gate | H8 | H13 | Requirement |
|---|---|---|---|
| G1 format | ✅ exact template, finite [0,1] inside, NaN outside | ✅ same | exact grid/CRS/transform/range |
| G2 uniqueness (343/343 corpus) | ✅ max \|r\| 0.0784, Jaccard 0.0386, containment 0.1057 | ✅ (re-audited) | \|r\|<0.85, J<0.50, C<0.80 |
| G3 blocked SGMC-proxy holdout | ❌ pooled DTI **0.003371** | ❌ pooled DTI **0.019677** | > 0.091550 (H33-B2), ≥2× random, > gradient |
| G4 repaired CAT-HID holdout | ✅ **0.023673** vs H33-B2 0.003590 | ✅ **0.021464** vs 0.003590 | strictly better |
| G5 LM-calibrated (diagnostic only) | 0.088394 | 0.053147 | not a gate (circular) |

Controls for G3: mass-matched random fields (seeds 40–42) reach 0.017663 (H8 budget) /
0.042564 (H13 budget); TMI-gradient top-K reaches 0.003495 / 0.007930. The circular SGMC-derived
incumbent scores 0.835907 on the same proxy (upper bound of SGMC-matching; disclosed, not a
target).

## Interpretation, without spin

1. **The Euler depth signal has real fault-finding skill** — both candidates beat H33-B2 by ~6×
   on the repaired catalogue-component holdout (G4). This is the first positive non-circular
   holdout result for the Euler family in this repository.
2. **…but that skill does not transfer to the off-catalogue SGMC proxy** (G3). The proxy's truth
   set is SGMC-derived; fields that do not imitate SGMC geometry top out far below the
   SGMC-lineage incumbent. H13's 5.8× improvement over H8 confirms the value-dilution diagnosis
   and shows the remaining gap is *placement*, not emission form.
3. **Random dots at equal mass beat both Euler fields on the proxy.** On the proxy, this says
   Euler crests are not preferentially located on off-catalogue SGMC traces — consistent with
   session 1's AUC 0.534 finding. On the live board, arrangement around true hidden traces
   (unknowable locally) is what the 0.24–0.28 family exploits; the proxy cannot rank that.
4. Two latent defects in the never-used CAT-HID instrument were found and repaired
   (hidden set included in its own exclusion flank; background mass in the hide fraction).
   Every historical CAT-HID number in this repository was a zero caused by these defects; none
   ever drove a decision.

## Why H33-B2 (0.2778) scored highest — measured answer

Byte-exact set comparison proves H33-B2 = H27-4 minus 2,545 pixels at ≤200 m from the catalogue.
Under `DTI = T/(0.2T+0.2FP+0.8G)` a removal that loses L of TP credit and saves S of FP mass helps
iff `(1−0.2D)·L < 0.2D·S`; at D = 0.2778 that is ≤0.0588 TP credit per FP unit (~17:1). The gain is
false-positive pruning of catalogue-flank mass, not new discovery. Whether 0.2778 can be beaten
locally is unknowable; the leaderboard does not identify TIFFs, and owner-reported file↔score
pairings remain unauthenticated. See
[`docs/data/h33-measured-analysis.json`](../data/h33-measured-analysis.json).

## Artifacts

| File | SHA-256 | Positive px | Mass Σp |
|---|---|---|---|
| [`gemsdoe40-h13-crest-binary-20261006-a5d5b80a8476.tif`](../downloads/gemsdoe40-h13-crest-binary-20261006-a5d5b80a8476.tif) | `c202a579ce33ea385a117f233857714150f8d68a5d6e8e46ef820b428e08ff33` | 21,041 | 21,041 |
| [`gemsdoe40-h8-tracelock-depthkde-20261006-373fa53b12e9.tif`](../downloads/gemsdoe40-h8-tracelock-depthkde-20261006-373fa53b12e9.tif) | `6c32147db39cd34e77baf31f6cb585165bca1d042cfa599d5d027cb140e0c195` | 213,616 | 7,993.2 |

Both are format-valid ([0,1] finite inside the 5,167,373-pixel footprint, NaN outside), unique
versus the full 343-artifact corpus, and HOLD under their own frozen gates. Receipts:
[h13-generation](../data/h13-generation.json), [h13-format](../data/h13-format.json),
[h13-uniqueness](../data/h13-uniqueness.json), [h13-validation](../data/h13-validation.json),
[h8-generation](../data/h8-generation.json), [h8-format](../data/h8-format.json),
[h8-uniqueness](../data/h8-uniqueness.json), [h8-validation](../data/h8-validation.json).

## Next experiment (registered, not implemented)

**H9 — blind basement-flexure faults** using bands never consumed by any arm
(`depth_to_base_surf` 15, `cond_surf` 17, `iso_grav_anom_hg` 18): a basement-depth flexure ridge
co-located with a conductivity break and gravity horizontal-gradient ridge marks buried normal
faulting under basin fill — the structure class a surface compilation least contains. Requires a
band-semantics audit first. Secondary: H10 geodetic strain-rate lineaments (bands 4/7/8, cheap),
H11 seismic-corridor × Euler intersection, H12 ComCat mechanisms (blocked by sandbox egress).
