# Candidate hypotheses H8–H12 and frozen H8 implementation — 2026-10-06 UTC (session 2)

**Registration status.** This slate was written and frozen **before** any H8 candidate code or
output raster existed. It is an append-only supplement to [`hypotheses.md`](hypotheses.md) and
[`preregistered-hypotheses-20261006.md`](preregistered-hypotheses-20261006.md). It does not revise
the H1/H2/H2-B/H4-A/H4/H7 results, and it does not re-open H4's HOLD.

## Decision target and evidence limits (unchanged)

Official task: predict faults that may indicate geothermal resources; metric is the 300 m
distance-weighted Tversky index, alpha = 0.2, beta = 0.8
([problem page](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/#performance-metric)).
Known USGS/INGENIOUS pixels are masked pixel-exactly when scoring
([staff clarification](https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516)),
and staff stated a "new fault" may include newly mapped geometry (continuations, splays, parallel
strands) of existing fault systems ([topic 11536](https://community.drivendata.org/t/11536)).
The only local truth-bearing instrument is the owner-derived SGMC proxy
(`26d142c4…`); its best comparator is itself SGMC-derived, so any proxy result is a **diagnostic**,
never an organizer score and never proof of hidden-set performance.

### Measured answer to "why did H33-2-B2 (0.2778) score highest?"

Byte-level set comparison in [`docs/data/h33-measured-analysis.json`](../data/h33-measured-analysis.json)
proves `h33-2-b2` is **exactly** the 40,199-pixel `h27-4` raster with the 2,545 pixels at ≤200 m
from the catalogue removed; it adds zero pixels. Under
`DTI = T/(0.2T + 0.2FP + 0.8G)` the marginal economics at DTI D say a removal that loses L of TP
credit and saves S of FP mass helps iff `(1 − 0.2D)·L < 0.2D·S`; at D = 0.2778 one may trade up to
0.0588 of TP credit per unit FP saved (~17:1). So the score gain is **false-positive pruning of
catalogue-flank mass**, not discovery. Whether any file can beat 0.2778 locally is unknowable: the
hidden labels are off-catalogue by construction, and the leaderboard does not identify TIFFs. What
*is* actionable: emit trace-continuous, depth-consistent, off-catalogue support and avoid diffuse
mass — exactly the two failure modes diagnosed in H4-A (KDE positive on 43.6 % of footprint,
AUC 0.534 vs catalogue) and in the H4 continuous field (proxy DTI 0.0069).

## Ranked hypotheses (registered before implementation)

### Rank 1 — H8: trace-locked Euler depth-consensus emission ("striking consensus")

- **Layers.** Hash-verified H4 Euler solution cloud derived from `tmi` (band 14, via
  `contact_euler.FIELDS` upward-continued 200 m) and `iso_grav_anom` (band 13, first vertical
  derivative as a top-edge approximation), SI = 0 with free contact offset A (Reid et al. 1990).
  Footprint from `sample_submission.tif`; exact catalogue pixels from `labels.tif`.
- **Physical signature.** A real near-surface fault produces *along-strike persistent*, shallow,
  mutually depth-consistent Euler solutions; noise and lithologic edges produce scattered or
  depth-incoherent ones. H8 therefore (a) keeps the H4 cluster-weighted KDE, then (b) restricts
  emission to orientation-NMS ridge crests of the combined field and (c) requires an oriented
  along-strike corridor (±800 m strike × ±300 m across) to contain ≥3 weighted solutions with
  depth standard deviation ≤ 445 m (normal-equivalent of MAD ≤ 300 m). Emission value is the
  combined consensus field times the depth-consensus factor, lightly smoothed (σ = 1 px) so the
  output stays continuous.
- **Why it can catch faults missing from USGS/INGENIOUS.** Continuation strands, splays and
  concealed extensions of fault systems keep coherent 3-D source geometry that single-window edge
  maps discard; exact catalogue pixels are zeroed, so nothing re-draws mapped traces. Staff's
  "new fault" definition explicitly admits such geometry.
- **Difference from anything in this repo or the 343-artifact corpus.** H1/H4/H7 emitted the raw
  KDE blob (H4: 865,145 positive cells, diffuse); H4-A showed that blob is not spatially
  selective. No repo arm and no corpus artifact applies orientation-NMS trace locking plus
  along-strike depth-consensus gating to an Euler cloud. The solver settings themselves are the
  frozen, already-validated H4 settings — novelty is claimed only for the emission geometry.
- **Expected DTI effect / cost.** Direction: converts the Euler family's proven depth signal into
  the trace-like arrangement that every ≥0.24 live artifact in the corpus uses. Magnitude unknown;
  no numeric forecast. Cost: low–medium (cloud already computed and committed; new code is
  vectorized numpy/scipy).

### Rank 2 — H9: blind basement-flexure faults (unused subsurface bands)

- **Layers.** `depth_to_base_surf` (15), `cond_surf` (17), `iso_grav_anom_hg` (18), `det_elev` (12).
- **Signature.** A flexure/step ridge in basement depth co-located with a conductivity break and a
  gravity horizontal-gradient ridge = buried normal fault under basin fill — the class of structure
  a surface compilation is least likely to contain.
- **Why off-catalogue.** Blind faults have no mapped surface trace by definition.
- **Difference.** No repo arm has ever consumed bands 15/17. This is the previously registered
  H4-E/H5 idea; registered here as the next experiment if H8 is evaluated.
- **Cost.** Medium; the layers are smooth and need a semantics audit (`cond_surf` units) before a
  step detector can be trusted.

### Rank 3 — H10: geodetic strain-rate lineaments

- **Layers.** `geod_shearrate` (7), `geod_dilaterate` (8), `geod_2ndinv` (4).
- **Signature.** Oriented structure-tensor lineament extraction on strain-rate fields; relay
  segments placed between, not on, catalogue traces.
- **Why off-catalogue.** Localized interseismic strain marks actively slipping strands that static
  compilations may not contain.
- **Difference.** Geodetic bands have entered prior arms only as ML features; no arm extracted
  oriented lineaments from them. Partially overlaps registered H6, but without the seismic-
  coincidence conjunction and with a simpler frozen emitter.
- **Cost.** Low. Risk: the grids may be too smooth to carry lineament-scale signal.

### Rank 4 — H11: seismicity-corridor × shallow-Euler intersection

- **Layers.** `deq_n100a15` (10), `ieq_n100a15` (16), H4/H8 Euler support.
- **Signature.** Cells simultaneously near recorded seismicity and under shallow depth-consistent
  Euler support mark active unmapped strands.
- **Difference.** Event layers were used as features before, never as an intersection gate with
  depth estimates.
- **Cost.** Low–medium; expected effect spatially limited to seismically active subregions.

### Rank 5 — H12: ComCat focal-mechanism orientation fields

- **External source needed.** USGS ComCat FDSN event web service
  `https://earthquake.usgs.gov/fdsnws/event/1/` (free, official).
- **Status.** `earthquake.usgs.gov` is outside this sandbox's egress; coverage gate (≥50 usable
  mechanisms over ≥8 blocks) unverifiable here. **Blocked — not viable in this session.**

**Ranking by expected DTI improvement per unit cost:** H8 > H9 > H10 > H11 > H12(blocked).
Only H8 is implemented this session.

## Frozen H8 implementation (all constants fixed before any scoring)

1. **Cloud input.** Read `docs/downloads/h4-euler-solutions.csv.gz`; fail closed unless SHA-256 is
   `6bed30b224981ef1064256d696c585b76f438011bbae6c0a16c00e6e6f5922ba` and 46,656 records parse.
   Easting/northing → fractional (row, col) with the exact sample transform. Retain records with
   `cluster_weight > 0` only. No re-solving; solver settings remain the frozen H4 `SETTINGS`.
2. **Per-family KDE.** Bilinear splat of weighted points per family; Gaussian σ = 2.0,
   truncate = 4.0; 99.5-percentile clip to [0,1] (frozen H4 `family_kde` semantics). Families:
   `tmi` (magnetic), `iso_grav_anom` (gravity). If either family has zero weighted points, stop.
3. **Consensus field.** `C = (M + G + sqrt(M·G))/3` (frozen H4 `combine_families`), zero outside
   the footprint.
4. **Orientation.** Structure tensor of C with integration σ = 3.0 px after inner σ = 1.0;
   strike θ = 0.5·atan2(2·Jxy, Jxx − Jyy) + π/2, quantized to 4 bins (0°, 45°, 90°, 135° strike).
   Cells with isotropy `(λ1 − λ2)/(λ1 + λ2) < 0.15` or C = 0 get no orientation.
5. **Ridge (orientation NMS).** Per quantized gradient direction (perpendicular to strike), a cell
   is a crest if `C ≥ forward-neighbour` and `C > backward-neighbour` (plateau-safe one-sided
   tie rule), or C equals both and is the scan-order first of an equal plateau is **not** required
   — the one-sided rule already thins plateaus to 1 px ridges. Union over bins.
6. **Along-strike depth consensus.** Oriented box kernels, per bin: half-length 8 px along strike,
   half-width 3 px across (built as explicit dense kernels ≤ 17×17). Four convolutions per bin on
   the solution splat grids: raw count `n`, cluster-weight mass `Σw`, weight-weighted depth
   moments `Σwz`, `Σwz²`; weighted mean `μ = Σwz/Σw` and corridor standard deviation
   `sd = sqrt(max(Σwz²/Σw − μ², 0))`. A crest cell passes iff `n ≥ 3` and `sd ≤ 445 m`.
   Consensus factor `γ = exp(−sd/445)`.
   *Pre-scoring implementation correction (same session, before any scoring): moments are
   cluster-weight-weighted rather than unweighted. Unweighted moments let incoherent, low-weight
   solutions inside the corridor veto coherent shallow clusters; the frozen H4 cluster weight is
   exactly the mutual-consistency evidence the corridor is meant to aggregate, so it must weight
   the moments. Raw-count gate unchanged.*
7. **Confidence.** `F0(crest) = C(crest) · γ(crest)`; `F = gaussian_filter(F0, σ = 1.0,
   truncate = 3)`; normalize by the in-footprint maximum to [0,1]; zero exact catalogue pixels
   (`labels == 1`); NaN outside the sample footprint; float32. No top-k cut, no proxy, no holdout,
   no prior prediction enters construction. If no crest passes, stop without relaxing gates.
8. **Determinism.** Fixed integer seeds where any randomness would be needed (none expected);
   float32 storage; canonical pixel SHA-256 recorded.

## Frozen validation and promotion gates (all decided before scoring)

- **G1 format.** `validate_candidate`: single band, float32, EPSG:32611, exact sample shape and
  transform, finite in [0,1] inside, all-NaN outside. Fail closed.
- **G2 uniqueness.** Full 343-blob corpus present in `data/prior`; `contact_audit.audit` passes
  both the H4 rule (|Pearson| < 0.85, top-37,654 Jaccard < 0.50, top containment < 0.80) and the
  historic H2-B rule against every exact-grid prior; no byte/canonical hash equality.
- **G3 blocked holdout (primary).** Frozen SGMC proxy `26d142c4…`, 4×6 blocks, 3-cell guard,
  exact-label zeroing (`score_array_on_proxy`): pooled DTI **strictly greater than 0.091550**
  (H33-B2 measured on this same instrument), and ≥ 2× the same-mass random control, and strictly
  above the same-mass gradient-top-K control.
- **G4 catalogue-component holdout (non-circular check).** CAT-HID mean strictly greater than
  H33-B2's CAT-HID mean, computed in this session with the same frozen `load_cat_hidden`.
- **G5 diagnostic only.** LM-calibrated mean reported; it is circular (incumbent and proxy are
  both SGMC-derived) and is explicitly **not** a gate.
- **Decision.** All of G1–G4 pass → status **SLOT-ELIGIBLE CANDIDATE** (owner decides; no upload
  happens from this repository). Any fail → **HOLD — research only**, named honestly. No
  post-score re-tuning of H8 constants; a variant would be a new preregistered hypothesis.

## What this session will *not* claim

No organizer score, no causal explanation of hidden labels, no "beats 0.3195" forecast. Proxy
passes are labelled proxy diagnostics. Owner-reported file↔score pairs remain unauthenticated
attributions under the standing integrity rule.

## Addendum — H13 registration and two instrument repairs (same session, 2026-10-06 UTC)

### Instrument defects found and repaired (before any H13 scoring)

While executing gate G4, `src/gems40/instrument.py::load_cat_hidden` was found defective in two
ways, both making its truth set wrong (the second made it exactly empty, which is why every
historical CAT-HID number in this repository is a zero):

1. `visible` included the hidden components themselves (`hidden ⊂ labels ⊂ visible`), so the
   3-px flank exclusion deleted the entire hidden set.
2. the 20 % hide threshold used `sizes.sum()` including the background class, selecting *all*
   catalogue components instead of ~20 % of catalogue pixels.

Both repairs restore the documented intent ("20 % of catalogue fault components are hidden and
scored in an eroded quadrant with a 1.5 km collar"). CAT-HID now yields 10,349 truth pixels across
four folds (17 % of catalogue after flank/domain exclusions). CAT-HID remains a diagnostic
instrument; it was never a decision instrument in this repository. All CAT-HID comparisons in this
session are re-measured under the repaired instrument for every candidate and baseline.

### H8 result recorded at the frozen gates (before H13 registration)

G1 format pass; G2 uniqueness pass (343/343 blobs; max |Pearson| 0.0784; max top-37,654 Jaccard
0.0386; max containment 0.1057; zero near-duplicates); **G3 fail** — pooled DTI 0.003371 vs
H33-B2 0.091550 and vs the mass-matched random control 0.017663. H8 is therefore **HOLD** under
its own frozen rules; no H8 constant is retuned.

### Rank 1.5 — H13: value-concentrated crest emission of the same depth consensus

**Metric-derived motivation (registered before any H13 scoring).** The official metric charges
`alpha * p(x)` per predicted pixel and grants at most `p(x) * k(d)` of credit per pixel, both
linear in the stored value. For any fixed value mass `Σp`, concentrating the mass into 1.0
predictions on the most fault-like cells maximizes attainable TP credit per unit FP charge;
max-normalizing a smoothed field instead dilutes mass into flank pixels worth at most their small
value. This is arithmetic from the metric definition
([problem page](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/#performance-metric)),
not a proxy observation.

- **Layers / physics.** Identical to H8 stages 1–6: frozen H4 cloud, per-family KDE, consensus
  field, structure-tensor orientation, orientation-NMS crests, along-strike depth consensus with
  the exact H8 constants. Nothing is re-derived or re-tuned.
- **Emission (the only difference).** Binary: 1.0 on crest cells that pass the depth-consensus
  gate, 0.0 elsewhere inside the footprint, exact known pixels zeroed, NaN outside. No smoothing,
  no dilation, no top-k budget: the crest-passing count is whatever the frozen physics gives.
- **Difference from corpus.** No artifact emits Euler depth-consensus crests as binary support;
  the nearest corpus members are topographic dotted ridges (different physics) and the H4/H8
  continuous KDEs (different emission).
- **Expected effect / cost.** Tests the single diagnosed weakness of H8 at zero new degrees of
  freedom. Cost: minutes; reuses all H8 code paths.

### Frozen gates for H13 (identical structure to H8)

G1 format; G2 full-corpus uniqueness (same thresholds); G3 pooled proxy DTI strictly greater than
0.091550, at least 2× the best mass-matched random control (mass = H13 crest count), and strictly
above the gradient-top-K control; G4 repaired CAT-HID mean strictly above H33-B2's repaired
CAT-HID mean; G5 LM diagnostic only. All pass → SLOT-ELIGIBLE CANDIDATE. Any fail → HOLD. No
post-score tuning; another variant would be a new preregistered hypothesis.
