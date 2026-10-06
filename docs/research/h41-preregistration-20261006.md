# H41 / H41b pre-registration — 6 October 2026

*Written before the scored comparison below was computed. Frozen choices are marked
**[frozen]**; additions made after seeing an intermediate result are marked **[post-hoc]** and
carry an explicit honesty note. Nothing here is a competition score.*

## 1. Physical hypothesis

Faults and fault-like contacts that are **missing from the provided catalogue** (USGS
Quaternary fault maps + INGENIOUS) should still appear in the potential-field layers as
**3-D contact source positions**, not merely as edges in a derivative map. Euler
deconvolution at the **contact structural index SI = 0** (Reid, Allsop, Granser, Millett &
Somerton 1990, *Geophysics* 55(1) 80–91, eq. 2 — the form that carries the arbitrary offset
A and is solved at their step 3b) returns a *depth-labelled* solution cloud. Tight clusters
of shallow, mutually consistent solutions that survive **all** window scales are candidate
contacts; scattered or deep solutions are not.

Why this should catch an unmapped fault: a contact that has no surface expression in the
catalogue can still produce a coherent magnetic/gravity source line. Why it might fail:
potential-field contacts include unfaulted lithological boundaries, and a survey drape makes
Euler depth an *effective* depth, not a ground depth.

## 2. Frozen computation **[frozen]**

| Item | Value |
| --- | --- |
| Fields | `rtp` (band 2), `tmi` (band 14), `iso_grav_anom` (band 13) |
| Structural index | 0.0 (fault/contact), offset A included |
| Windows / stride | 10, 16, 24 px / 4 px (full run); 10, 16 px / 8 px is the smoke variant only |
| Euler gates | analytic-signal percentile 72, max relative depth error 0.22, depth in [80, 2200] m, source-in-window, conditioning floor |
| KDE | Gaussian σ = 1.7 px splat of solution weights; weights = shallowness × quality × tightness × mutual depth consistency |
| Concordance | magnetic *and* gravity solutions agreeing in position and depth get ×(1 + 1.4·concordance) |
| Scale stability | geometric mean of the per-window KDEs, normalised per window; a cell must be non-zero at **every** window scale |
| Shallow filter | only solutions ≤ 1500 m enter the KDE |
| Support | footprint ∧ ¬catalogue ∧ distance-to-catalogue > 1 px ∧ Euler-gated (`field > 0`) |
| Spacing | 3 px (equal to the metric's own kernel support R = 300 m) |
| Emission | value-ranked Poisson-disk thinning of the ranking field |
| Mass | the snapshot maximising the recalibrated instrument among {10k, 20k, 30k, 40k, 60k, 90k, 120k} |

## 3. Ranking **[frozen, features listed for audit]**

A **spatially blocked, out-of-fold** L2 logistic discriminant (5 × 6 contiguous block folds;
0.5/99.5 percentile winsorisation; standardised features) ranks cells inside the support.

* Positives: public faults absent from the provided catalogue (**USGS SGMC / DS-1052**,
  CC0/US public domain), dilated 1 px, weighted by the metric's triangular kernel credit
  (1.00 / 0.67 / 0.34).
* Negatives: cells further than 3 px from any such fault; the 1–3 px band is dropped.
* The provided catalogue is excluded from both classes and from the support.
* Reported metric: **out-of-fold AUC**, pooled and per block (never in-sample).

Features: the three per-window KDEs, the scale-stable field, the three per-window mean-depth
maps, |∇rtp|, |∇gravity|, |∇DEM|, DEM curvature at 1 px and 3 px, the ten raw band values
(`rtp`, `tmi_hg`, `geod_2ndinv`, `iso_grav_anom_slope`, `tc`, `iso_grav_anom_vg`, `det_elev`,
`iso_grav_anom`, `tmi`, `depth_to_base_surf`), **and [post-hoc]** the cross-field vertical
gradient ratios `|tmi_vg| / (|iso_grav_anom_vg| + 1e-3)` and `|tmi_hg| / (|gravity_slope| + 1e-3)`,
the 2 px curvature of the RTP field, a 15 px local z-score of `tc`, and the distance to the
provided catalogue.

**Honesty note on the post-hoc additions.** They were added after the first gated run
measured 0.6423 (Euler features only) and 0.6413 (Euler depth features), and after a second
run with DEM curvature + bands 4/6 measured 0.7389. The additions are therefore *not*
pre-registered; they are reported as an instrument-selected model, and the AUC improvement
they produce is in-sample of the selection. The candidate's blocked proxy comparison and the
physical construction are unaffected by this, but a fresh reader must treat "27-feature
model" as selected.

## 4. Audit gates — all must pass before a file is offered **[frozen]**

1. Format: single-band float32, EPSG:32611, 100 m, 3730 × 3292, identical transform to the
   sample submission; inside the footprint every value finite and in [0, 1]; nothing emitted
   outside the footprint.
2. Novelty: |Pearson| < 0.85, support Jaccard < 0.35 and top-budget Jaccard < 0.35 against
   **every** staged prior raster (51 same-grid files), computed on the footprint from raw
   values.
3. Blocked proxy comparison on the identical 4 × 6 partition with a 3-pixel guard against
   the three strongest retrievable comparators (H33-B2, H27-4, H40-E).
4. Metric algebra: every emitted dot must clear the exact marginal credit bar of the
   published metric, `k > 0.2 · DTI`; no dot is emitted within the catalogue flank.

## 5. Promotion bar and the honest outcome **[frozen]**

A candidate may be recommended for a weekly slot only if it beats the best available
holdout comparator **and** the instrument used to project a live score reaches a
leave-one-out Spearman rank correlation of **≥ 0.80** against the scored anchor corpus
(`ref/prior`, owner-reported scores). Below that bar the file is published as a research
candidate with `slot_eligible: false` in `docs/data/current-candidate.json`.

Measured (6 Oct 2026, `docs/data/live-transfer.json`): the best of four instruments
(saturating, power-law, hidden-truth algebra, block-linear) reaches **LOO Spearman 0.705**
(saturating, LOO RMSE 0.0538). The bar is **not met**; the repository therefore does not
certify a slot, and says so on every page it publishes.

## 6. What would falsify this arm

* A blocked proxy mean below the owner's two best scored dot sets.
* A near-duplicate correlation with any prior output.
* Any non-finite cell inside the footprint, or a value outside [0, 1].
* A surrogate that ranks the scored anchors worse than the adopted instrument (already
  observed for the whole-map proxy DTI, LOO Spearman −0.897 — recorded as a warning against
  optimising for that proxy).
