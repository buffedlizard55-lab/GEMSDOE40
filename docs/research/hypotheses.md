# Preregistered research hypotheses — GEMSDOE40

**Lock date:** 2026-10-05 (UTC)
**Status:** Registered before any candidate-generation code or candidate TIFF was written.
**Decision principle:** *Maximize P(Win)*; *Own the Outcome*. No leaderboard score, holdout result, or submission is assumed.

This is a research register, not a claim that any method will win. The challenge target is fault probability, not geothermal-vent location. A potential-field contact is not automatically a fault, and a fault is not automatically a viable geothermal reservoir.

## Evidence and starting point

- The official task is to predict geologic faults that may indicate geothermal resources. The official metric is a continuous, 300 m triangular distance-weighted Tversky index with `alpha=0.2`, `beta=0.8`; see the [official problem and metric page](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/).
- The official rules require one single-band float32 GeoTIFF at 100 m in EPSG:32611, same bounds as the training data, with null/NaN outside and values in `[0,1]` ([September 2026 rules](https://docs.nlr.gov/docs/fy26osti/96647.pdf)).
- DrivenData staff clarified that the known USGS/INGENIOUS mask is **pixel-exact**, not a 300 m exclusion buffer; predictions adjacent to known traces are still penalized ([staff clarification](https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516/2)).
- Euler deconvolution is prior art in [GEMSDOE28](https://buffedlizard55-lab.github.io/GEMSDOE28/): its H31-1 reports TMI-only SI=0 solutions and depth clustering; H38-1 added heat-flow evidence. GEMSDOE40 therefore does **not** claim novelty for Euler itself. Any novelty claim is limited to the specifically preregistered joint-field/depth-consensus construction and must pass the raw-raster audit.
- The available off-catalogue comparison raster is an owner-derived mirror of an SGMC-derived product, not hidden challenge truth and not an independently rebuilt USGS product. It is a **proxy instrument only**.

## Ranked hypotheses (registered before implementation)

### Rank 1 — H1: joint magnetic–gravity Euler depth-consensus KDE

**Hypothesis.** Locations supported by mutually close, depth-consistent 3-D Euler solutions from the magnetic field and the vertical derivative of the gravity field are more fault-like than locations supported by either field alone. A local 2-D KDE of those matched solution points should improve precision at a fixed emission budget relative to the best prior output on the frozen off-catalogue proxy blocks.

- **Layers / data:** challenge `tmi` (band name, not assumed band number), `tmi_vg`, and `iso_grav_anom`; 100 m grid. Horizontal derivatives will be computed from the named scalar fields with the raster transform. The gravity vertical derivative will be derived from `iso_grav_anom` by a padded Fourier-domain potential-field derivative; the provided `iso_grav_anom_vg` band will not be trusted as a vertical derivative unless independently validated. The magnetic `tmi_vg` band will be used only if its Pearson correlation with a reflect-padded (128-cell) Fourier `+|k|` downward derivative of `tmi` is at least 0.95 and its least-squares scale ratio lies in `[0.8, 1.2]`; otherwise the derived band will be used. Horizontal derivatives use central finite differences at 100 m. For gravity, derive `Gz_down = F⁻¹(+|k|F(G))` and the Euler vertical derivative of `Gz_down` with the same operator; never substitute the supplied `iso_grav_anom_vg` unless a separate validation passes.
- **Physical signal:** a fault-zone density/magnetization boundary may produce related magnetic and gravity sources; Euler supplies location and depth hypotheses rather than only a 2-D edge. Cross-field agreement is evidence, not proof, because lithologic contacts and intrusions can create the same signals.
- **Frozen method:** 10×10-cell windows, 5-cell stride; SI=0 for magnetic contacts; SI=0 on first vertical gravity derivative as an approximation to a finite gravity contact (raw finite-contact SI is −1 in the cited idealization). Reject only on predeclared algebraic/physical checks: relative Euler fit residual >0.20, normal-matrix condition number >1e5, source outside the horizontal half-window, or depth outside 100–10,000 m. For each magnetic solution, consider gravity solutions within 300 m XY and 1,000 m depth; choose the gravity solution minimizing `(dXY/300)^2 + (dZ/1000)^2` (gravity solutions may be reused). Pair weight is `exp[-0.5((rm/0.20)^2 + (rg/0.20)^2 + (dXY/300)^2 + (dZ/1000)^2)]`. Deposit each pair at its XY midpoint, sum weights per cell, then apply a zero-padded Gaussian filter with sigma=2 pixels and truncate=3.0. Emit the top **45,962** eligible cells at binary 1.0 (the previously published H40-F emission count), zero elsewhere inside the template footprint, and zero on the exact known-label pixels; break exact ties by ascending row-major cell index. Fail rather than backfill if fewer than 45,962 eligible cells have positive KDE. No holdout labels, SGMC pixels, or leaderboard score enter feature construction or ranking.
- **Why it may find off-catalogue faults:** two potential-field physics channels constrain a 3-D source more strongly than a single TMI cluster, and the exact-label mask means the target is not simply to redraw the supplied traces. New fault geometry can still lie near a known trace; no distance buffer is applied.
- **Novelty boundary:** different from GEMSDOE28's TMI-only Euler clouds / heat-flow overlay because the decision variable is *cross-field, depth-consistent solution-pair density*, not a TMI point cloud or a generic edge threshold. This is a narrow method-level novelty claim, pending independent code review and a raw-output hash/correlation audit against every available prior submission raster.
- **Expected improvement:** directionally positive versus single-field Euler if the added gravity constraint removes magnetic-only artifacts; magnitude is unknown and no numeric score is forecast. The likely failure modes are mismatched source physics, derivative noise, and gravity/magnetic responses to non-fault lithologic contacts.
- **Cost:** medium–high; two full-grid derivative pipelines, windowed least-squares inversions, spatial pairing, KDE, and full-raster audits. No new external data required.

### Rank 2 — H2: upward-continuation persistence of magnetic Euler solutions

**Hypothesis.** Source solutions that persist in location and depth after successive upward continuation of `tmi` are more likely to represent coherent deeper structures than shallow noise or processing artefacts. A persistence-weighted solution cloud should improve over single-height TMI Euler clustering.

- **Layers / data:** named challenge layer `tmi` plus its mathematically derived horizontal and vertical derivatives; upward continuations fixed at 0, 500, 1,000, and 2,000 m.
- **Physical signal:** upward continuation attenuates short-wavelength components while preserving the geometry of deeper sources; a real source should not migrate arbitrarily as height changes.
- **Why it may find off-catalogue faults:** it tests for coherent subsurface structure without copying the known-fault mask and can suppress shallow cultural or near-surface magnetic noise.
- **Novelty boundary:** GEMSDOE28 documents a single-height Euler solution cloud; the proposed persistence trajectory across four heights is the distinct element. It remains in the Euler family and would not be described as wholly novel Euler theory.
- **Expected improvement:** modest/uncertain; likely better stability and lower false-positive density, with a risk of erasing shallow fault sources. No numeric improvement is claimed.
- **Cost:** high; four FFT continuations and four repeated Euler solves, plus depth/position trajectory matching.

### Rank 3 — H3: focal-mechanism-constrained fault-plane orientation

**Hypothesis.** USGS ComCat focal mechanisms provide strike/dip/rake constraints that can distinguish active fault planes from unrelated geophysical lineaments. A fault-plane orientation likelihood, intersected with independent potential-field support, may improve off-catalogue fault ranking beyond the supplied earthquake-density layer.

- **Layers / data:** official USGS ComCat earthquake locations and focal-mechanism/nodal-plane fields, with challenge `deq_n100a15`, `geod_shearrate`, and potential-field lineaments as context. The API endpoint was reachable and returned 14,665 events for the preregistered Great Basin bbox/date/magnitude query; the number and spatial coverage of actual focal mechanisms have **not** yet been verified.
- **Physical signal:** focal mechanisms constrain the orientation of slip planes for recorded earthquakes; they do not directly map inactive faults or prove geothermal permeability.
- **Why it may find off-catalogue faults:** event mechanisms are independent evidence of fault orientation not represented by the smooth earthquake-density feature alone.
- **Novelty boundary:** the supplied earthquake-density band is a scalar activity proxy; the distinct signal would be mechanism-derived plane geometry. No claim of effectiveness until an official-source payload and coverage check are complete.
- **Expected improvement:** potentially useful only for the seismically active subset; broad-region effect is unknown and could be negligible.
- **Cost / feasibility:** high and currently **blocked from viability**. The ComCat API is obtainable, but the mechanism count/coverage has not been checked; this hypothesis must not be implemented or called viable until that check passes. No additional service or paid dataset is proposed.

## Frozen validation and promotion rules

1. **Only H1 is the first implementation candidate.** H2 and H3 will not be used to tune H1. H1's window, SI, fit, depth, pairing, KDE, and 45,962-cell emission settings above are locked before the proxy score is read.
2. **No label leakage.** Candidate construction uses the feature bands and sample footprint only. The known training labels are used only for exact-pixel suppression and to remove exact known pixels from the SGMC proxy truth; there is no 300 m buffer. SGMC is never used for fitting, pairing, thresholds, or emission ranking.
3. **Spatial holdout:** partition the common template into 4×6 contiguous blocks. Score only block interiors after a three-cell guard at each block boundary; use the official 300 m DTI definition, pooled by summing weighted TP/FP/FN across those held-out interiors. Report every block and the pooled result.
4. **Incumbent:** score every available, same-grid prior prediction TIFF from the public sibling-repository inventory on the exact same frozen blocks; select the highest pooled proxy DTI as the holdout incumbent. The SGMC raster is a derived proxy, so these numbers are not organizer scores and do not predict the private leaderboard.
5. **Promotion gate:** H1 is eligible for a weekly submission slot only if (a) its pooled proxy DTI strictly exceeds the incumbent, (b) it is higher in at least 18 of 24 blocks, and (c) no block is lower by more than 0.005 absolute. A failure is a valid negative result: retain the TIFF as research-only and explicitly say **DO NOT SUBMIT**. No weekly submission slot is spent otherwise.
6. **Independent artifact checks:** exact sample transform/CRS/shape/bounds/nodata semantics; single-band float32; finite `[0,1]` in-footprint values; null only outside; raw SHA-256 and canonical-pixel hash; Pearson correlation, binary-support Jaccard, and top-budget overlap against every comparable prior raster. A hash match or near-duplicate under the preregistered audit thresholds blocks any uniqueness claim.
7. **Near-duplicate rule:** uniqueness fails if any prior same-grid raster has the same raw SHA-256, the same canonical in-footprint pixel array, `|Pearson r| ≥ 0.995` over the sample footprint, nonzero-support Jaccard `≥ 0.90`, or top-budget support containment `≥ 0.90` (top budget is 45,962 positive cells). The inventory is deduplicated by Git blob but retains every observed path/ref; every comparable raster is checked.
8. **No outcome inflation:** public leaderboard score, proxy DTI, owner-reported score, and projected/modelled score are separate evidence classes. Only a row on the official leaderboard establishes an account-level public score; it does not identify a TIFF unless there is a direct submission receipt.

## Append-only amendment — operational lock for H2 after H1 feasibility stop

**Timestamp:** 2026-10-05 UTC, before H2 code and before any proxy holdout score was read. The first H1 implementation returned 2,305 accepted magnetic solutions and **zero** gravity solutions under the locked quality gates, so its pre-registered 45,962-cell KDE could not be formed. H1 is stopped; its thresholds will not be relaxed and it will not be backfilled. This is a feasibility result, not a holdout score.

If H2 is evaluated, its exact implementation is now locked as follows:

- Use named `tmi` only. Upward-continue the filled TMI grid at exactly 0, 500, 1,000, and 2,000 m using a 128-cell reflect pad and `exp(-|k|h)` in the Fourier domain; derive horizontal derivatives by 100 m central differences and downward-positive vertical derivative by `+|k|`.
- At every height use 10×10 windows, 5-cell stride, SI=0, condition number ≤1e5, relative Euler residual ≤0.20, horizontal source within ±500 m of its window centre, and source depth corrected to the ground plane as `d_ground = d_solution_from_observation − h`, retaining 100–10,000 m only. Use the same one-cell data-support halo as H1.
- Anchor matches on each accepted height-0 Euler solution. At each other height choose the solution minimizing `(dXY/300)^2+(dZ/1000)^2` among those with `dXY≤300 m` and `|dZ|≤1,000 m`; solutions may be reused. Require matches at at least 3 of the 4 heights (including height 0). Weight is `(n_hits/4) × exp[-0.5 × mean(residual²/0.20² + normalized_match_distance²)]`, with the height-0 match terms set to zero.
- Deposit at the height-0 source location, apply the same zero-padded Gaussian KDE (sigma 2 px, truncate 3), and emit the top 45,962 positive eligible cells with row-major tie breaking; exact known-label pixels are suppressed, no distance buffer is used, and failure to reach the fixed budget aborts H2.
- H2 inherits the same frozen 4×6 spatial-block proxy holdout, incumbent set, promotion gate, format validator, uniqueness thresholds, and **DO NOT SUBMIT** failure action. No H2 settings may be selected using proxy DTI.

This is an append-only contingency, not a rewrite of the original ranking. H3 remains blocked until its ComCat focal-mechanism coverage is measured.

### Append-only H3 source-coverage audit lock

Before querying event records, measure whether the official USGS ComCat source has adequate coverage for H3 using the sample footprint's geographic bounding box transformed from EPSG:32611: longitude −120.037171 to −116.140922, latitude 37.331193 to 40.727879. Query event origins from 1970-01-01 through 2026-10-05 UTC with no magnitude cutoff and `producttype=focal-mechanism`; then inspect returned event-detail products for at least one finite nodal plane strike/dip per usable mechanism. Do not use SGMC, the owner proxy, challenge fault labels, or leaderboard values in this source audit. H3 is considered source-coverage-feasible only if at least 50 usable mechanisms are distributed across at least 8 of the registered 24 spatial blocks; otherwise stop H3 without implementation. This threshold is a broad-region coverage screen, not a claim that earthquakes map inactive faults.

### Append-only H2-B: support-limited H2 output (locked before proxy scoring)

The strict H2 top-45,962 implementation has already demonstrated that only 17,425 valid positive-KDE cells exist. It aborts rather than adding unsupported cells, as preregistered. Before any holdout score is calculated, a separate **H2-B** output rule is registered: keep the same four-height Euler persistence field and KDE, but emit **all and only** in-footprint, non-known-label cells with KDE `> 0` at 1.0; emit zero at every other in-footprint cell; never add a top-K threshold or backfill. This is a deterministic support-limited scientific candidate, not a reinterpretation of H2's failed fixed-budget gate. It inherits the same format, uniqueness, 24-block holdout, and promotion checks. No holdout-label or SGMC value enters its support or ranking. If it does not beat the frozen incumbent under the inherited gate, it is marked **HOLD; DO NOT SUBMIT**.

### Append-only H3 source audit outcome — still blocked

The official ComCat count endpoint returned 34,716 event records with a `focal-mechanism` product in the locked footprint/time query, and 238 with magnitude ≥4.0. These product-filter counts use the newly frozen sample-extent query and are not directly comparable to the earlier 14,665-event availability query recorded in H3's original preregistration. One event-detail payload (`nc10085763`, M6.4) exposed two finite nodal planes. This confirms source access and that mechanism fields exist, but it does **not** verify the number of usable solutions across all events or their distribution over the registered 24 blocks. Direct sandbox Python HTTPS requests to the official host failed with TLS/SSL connection closure; the research web tool could read the count and sample event detail. The registered threshold of ≥50 usable mechanisms in ≥8 blocks therefore remains unverified; H3 is **not implemented and not called viable**. See [`h3_source_audit.json`](../data/h3_source_audit.json).

## Status changes

This register is append-only. After a result is available, append the result and any deviation in `docs/research/validation.md`; do not silently edit the registered thresholds or call proxy performance an official score.

## Append-only register — H40 session (2026-10-06 UTC)

Registered before the H40 emission rule was implemented. The brief demands three to five
candidate geological hypotheses, ranked by expected gain against cost, each naming its layers,
its physical signature, why it could mark a fault that the USGS/INGENIOUS catalogue misses, and
how it differs from everything already in this repository.

| rank | hypothesis | layers | physics | why it can catch an unmapped fault | novelty in this repo | cost | outcome |
|---|---|---|---|---|---|---|---|
| 1 | **H40-A Euler SI=0 depth-cluster crest** | `rtp` (B2), `tmi` (B14), `iso_grav_anom` (B13) | Reid et al. (1990) deconvolution at the structural index of a *contact* (SI = 0): a fault offsets magnetisation and density, so its solutions cluster tightly in (x, y, depth) and agree in depth with each other | the catalogue is compiled from field mapping and literature, never from gridded potential-field deconvolution; a buried contact whose surface trace is unmapped still produces a solution cluster | maximum support-Jaccard 0.02 and max \|Pearson\| 0.06 against all 26 staged priors; no gradient/ridge/curvature field in the corpus is built this way | 12 min | **implemented — does not pay (see below)** |
| 2 | H40-B magnetic × gravity Euler depth conjunction | `rtp`, `tmi`, `iso_grav_anom` | keep only magnetic solutions that have a gravity solution within 3 px **and** within 35 % relative depth: a fault that offsets susceptibility *and* density | gravity Euler solutions are the scarce measurement (10,255–34,450 accepted here vs 40,130–54,060 magnetic); a conjunction is far more selective than either family | no staged prior uses gravity Euler deconvolution at all | 1–2 sessions | not implemented — cheapest next test |
| 3 | H40-C strike-continuation on the catalogue flanks | any structural layer | expert-added faults are usually along-strike extensions of mapped traces | — | — | 1 session | **refuted** by the paired live measurement: deleting exactly the 6,436 dots at 100–200 m from mapped traces raised the live score 0.2600 → 0.2778 (+0.0178), i.e. those dots earn ≈ 0 credit |
| 4 | H40-D depth-sliced emission | same as H40-A | emit only solutions shallower than 400 m; the sibling repository measured off-catalogue SGMC enrichment of 4.4× at 50–100 m falling to ~1.0× below 800 m | shallow contacts are the ones expressed as surficial faults | depth weighting is already inside `solution_weights`; a hard slice is new | hours | tested only as a soft gate (250 m, 120 m): w 0.0520 → 0.0536, predicted 0.1575 → 0.1609 — within noise |
| 5 | H40-E hidden-truth density inversion from the public score record | the 12 organizer-scored artifacts | treat the live scores as constraints on a blockwise hidden-truth density and weight the emission by the posterior | the score record is the only direct measurement of where the hidden truth actually is | nothing in the corpus uses the score record as a spatial instrument | 1 session | not implemented; needs multiple-comparison control (12 constraints) |

### H40 outcome (measured, not asserted)

- **Euler SI = 0 runs on both families with the locked gates** (stride 4, analytic-signal
  percentile 72, max relative depth SE 0.22, depths 80–2200 m, ~8 s per configuration):
  `rtp` 54,060 / 47,174 / 40,130 accepted solutions at windows 8 / 12 / 16 (median depth
  194 / 192 / 189 m); `iso_grav_anom` 10,255 / 23,339 / 34,450 (median depth 355 / 534 / 669 m).
  This is a *different* configuration from the preregistered H1 stop (which had zero accepted
  gravity solutions under H1's own gates) and does not overturn that record.
- **The emission rule was selected by the only live-anchored instrument available.** The
  saturating model fitted by the sibling repository to 12 organizer-scored artifacts is
  reproduced here to 5 × 10⁻⁵ on its input quantity `w` (mean per-dot kernel credit against the
  off-catalogue surrogate truth) and to 2 × 10⁻⁴ on its published predictions (9 anchors present
  locally). Scanning 5 spacings × 7 budgets chose **spacing 2.5 px, mass 60,000**, with
  `w = 0.0544` and predicted live score **0.1625** (break-even marginal `w` 0.0537 — the chosen
  emission sits exactly at its own optimum).
- **The pre-existing promotion gate is invalid and is not used to promote this artifact.** Its
  comparator is the highest pooled SGMC proxy score among format-eligible priors, which is the
  *blind spacing-5 lattice* (proxy 0.2679 pooled / 0.1650 mean, owner-reported live **0.0904**),
  while the incumbent that leads the live record (0.2778) scores 0.1096 / 0.0722 on the same
  instrument. A gate whose best comparator is the worst live artifact cannot rank candidates;
  the evidence is recorded in `docs/data/gate-defect-20261006.json` rather than silently dropped.
- **Verdict: HOLD — DO NOT SUBMIT.** Two live-anchored instruments disagree about the magnitude
  (saturating instrument 0.163; the 4-quadrant prevalence-calibrated LM instrument 0.202,
  which is 0.066 below the same instrument's incumbent). Neither reaches the incumbent's
  reported 0.2778, so no submission slot is justified, and the artifact ships as a research
  record only.
- **Corrected calibration of the metric record.** With the corrected `w` values, the two files
  that bracket the catalogue-flank pruning give, jointly, hidden-truth credit
  `T ≈ 5,220` at `M = 44,090` for the 0.2600 artifact (surrogate-to-live credit ratio
  `c ≈ 0.54`) and `0.8·N + 0.2·(T − C) ≈ 11,261`, i.e. `N ≈ 12,800–14,100`. An earlier
  session's claim that adding Euler dots to the incumbent "pays +0.008 to +0.041" used the
  *unweighted* surrogate credit, which over-rewards breadth (the blind lattice earns 0.142
  surrogate credit per pixel against 0.0505 per dot on the instrument definition). At the
  corrected bar, Euler dots at `w = 0.053` do **not** pay (the incumbent's own operating point
  requires marginal `w > 0.099`).

### Addendum 2026-10-06 (later in the same session): instrument audit and H40-B

**Instrument audit — no available instrument can rank the scored artifacts.**
The 16 organizer-scored artifacts that exist as rasters in this corpus were each
scored with the two instruments the site had been using, plus the raw
off-catalogue SGMC credit density `w` (full table:
`docs/data/instrument-audit-20261006.json`):

| quantity | Spearman vs the 16 organizer scores |
|---|---:|
| 4-quadrant LM instrument (`lm_calibrated`) | **+0.10** |
| emitted mass `M` | **−0.676** |
| off-catalogue SGMC credit per dot `w` | **+0.676** |

Two decisive counterexamples destroy the instruments as ranking tools:

1. **Mass-matched control pair.** `p34-scatter-q50` and the incumbent
   `gemsdoe32-h33-h33-2-b2` both emit exactly 37,654 dots. The control holds
   `w = 0.1896` of SGMC credit per dot, the field holds `w = 0.1776` — 7 % *less* —
   yet the organizer scores are 0.0778 and 0.2778. A 3.6× inversion at equal mass
   means SGMC credit density carries no information about placement quality.
2. **Blind lattice.** `p13-lattice-s5-v2` is a blind 5-px lattice with no
   geological input. The LM instrument ranks it **first of all 16** (0.4366) while
   its organizer score is third from the bottom (0.0904). The LM instrument's
   +0.10 rank correlation is the whole story: it rewards breadth, because its
   surrogate truth set (62–80 k px) is 5–6× denser than the hidden truth
   (`N ≈ 9–14 k px`).

The earlier statements in this register and in the artifact receipt that quoted
"predicted live 0.1625" and "LM 0.2015 vs incumbent 0.2679" are therefore
**withdrawn as promotion evidence**: the saturating model's input `w` is defined
on a different surrogate than the one it was fitted with, and refitting the same
functional form with the locally defined `w` gives leave-one-out Spearman 0.47
with the two best artifacts inverted. The H40 artifact stays **HOLD**, now for a
stated and measured reason rather than an instrument number.

**H40-B — magnetic × gravity Euler depth conjunction (measured, not refuted).**
Requiring every magnetic SI = 0 solution to have a gravity solution within 3 px at
±60 % relative depth keeps 11,848 of 285,536 magnetic solutions (4.1 %; 2,163 /
2,047 / 1,654 for `rtp` at windows 8/12/16 and 2,221 / 2,026 / 1,737 for `tmi`).
The resulting field, thinned at 30,000 dots (spacing 3.0), reaches
`w_offcat = 0.1333` against `0.1277` for the all-solution field at the same mass
(+4 %), and 0.1508 at 15,000 dots — the highest credit density measured in this
field family. Under the audit above this is *not* promotion evidence, but it is
the cheapest next measurement: H40-B is retained as the top-ranked untried
hypothesis, and the next session's mandatory step is a **truth model validated by
leave-one-out ranking of the 16 scored artifacts (target Spearman ≥ 0.8)** before
any emission decision.
---

## Append-only amendment — H4 slate, registered 2026-10-06 (UTC), before any H4 code or candidate raster was written

**Why a new slate.** The frozen all-prior promotion gate cannot be passed by an honest independent candidate: the frozen incumbent (`gapfinder-v2-sgmc-gap`, pooled proxy DTI 0.835907) was constructed from the SGMC-derived proxy pixels themselves, so scoring against that proxy is circular. That limitation is documented in the README and the H2-B report. The H4 slate therefore (a) keeps the frozen gate as a *diagnostic that is expected to fail*, and (b) registers one new, non-circular comparative instrument that can rank *our* candidate against the whole prior corpus on the same blocks with the same rule. Both results are reported side by side; no threshold in the frozen gate is relaxed.

**New comparative instrument (registered now, before use).** "Catalogue-block instrument": the organizers' own existing-fault raster (`data/labels.tif`, 60,988 positive 100 m cells inside the sample footprint) is the truth; the scored domain is the sample footprint; the holdout is 4 rows × 6 columns of contiguous blocks with a 3-cell guard on every internal edge (identical geometry to the frozen instrument); the metric is the published distance-weighted Tversky index (alpha=0.2, beta=0.8, R=3 px). Every prior GeoTIFF that shares the sample grid is scored by exactly the same code. **Known limitation, registered up front:** the hidden target is *off-catalogue* by construction, so this instrument measures fault-finding skill against a published, partially misaligned catalogue, not hidden-set performance. It is a comparator, not a forecast.

### H4-A (rank 1 — implemented this session): Euler SI=0 depth-cluster contact network with metric-optimal dotted emission

- **Layers / data.** Challenge bands `rtp` (2) and `iso_grav_anom` (13); the gravity vertical derivative is derived in the Fourier domain (`+|k|` operator, Blakely 1995 eq. 12-8), never taken from the supplied derivative band. Surface-expression concurrence uses the 3DEP 1 m DEM–derived 100 m scarp-feature stack (`external/dem/lidar_scarp_features_u8.tif`, 12 channels, USGS 3DEP, no use restrictions) plus challenge band `det_elev_slope` (19). Mapping-gap weighting uses the official USGS Quaternary Fault and Fold Database raster (`Qfaults_GIS.zip`, split by mapped scale).
- **Physical signature.** A *source locator*, not a threshold: sliding-window least-squares Euler deconvolution at structural index N=0 (fault-like contact; Reid, Allsop, Granser, Millett & Somerton, *Geophysics* 55(1):80–91, 1990) yields depth-labelled points; the candidate field is a kernel-density estimate of solutions that are (i) shallow, (ii) well fitted, (iii) locally tight in 3-D and (iv) mutually consistent in depth, with magnetic and gravity solutions that agree in (x, y, z) boosted. A multi-scale Hessian lineament response then converts the point-density field into *trace* geometry, because the target is a line raster, not a point cloud.
- **Why it should catch faults missing from the catalogue.** The catalogue is the USGS Quaternary fault compilation (this run measures 58,198 of the 60,988 catalogue cells inside the footprint to lie on USGS QFFD coarse-scale traces, 95.4 %) — i.e. it is a *coarse-scale, surface-based* compilation. (a) A potential-field contact locator images structural discontinuities including buried and intra-basin structures that a surface compilation does not carry; (b) 3DEP 1 m lidar is exactly the evidence class used for detailed fault mapping, so scarp-expressed traces that the 1:250k compilation omits are the expected gap; (c) requiring magnetic/gravity depth agreement and surface expression suppresses the lithological contacts that dominate either channel alone.
- **How it differs from everything already implemented in this repository.** H1 (joint mag–gravity pairing) stopped with zero accepted gravity solutions under its locked gates and never produced a raster; H2/H2-B deconvolved `tmi` only, at four upward-continuation heights, with a fixed top-budget or support-limited emitter and no surface-expression term; `scripts/run_euler_depth_cluster.py` produced a *continuous KDE* whose positive support was every cell above a floor quantile, with no trace extraction, no dot decimation and no terrain concurrence. H4-A adds: multi-window (8/12/16 px) SI=0 clouds from both fields with a 3-D depth-consistency *gate* (not only a weight), a multi-scale Hessian lineament transform of the Euler field, a 3DEP/DEM concurrence term, mapping-gap weighting, and an emission rule derived from the exact metric identity `DTI = TPw / (0.2*(TPw+FPw) + 0.8*|G|)` with non-maximum-suppressed dots at the registered 3 px spacing.
- **Pre-registered emission rule.** Dots are placed by greedy non-maximum suppression on the final field with a hard minimum separation of 3 px (300 m), the spacing family that the best live-scored artifact used (median nearest-neighbour distance 3.0 px measured on `gemsdoe32-h33-h33-2-b2-…-nan.tif`, 37,654 isolated 1-px components) and that `docs/data/emission_theory.json` measures against the catalogue. Emission is suppressed at every exact catalogue pixel and within 2 px of it (the prune that the only live-scored ablation in the family supports), and the number of dots is capped at the registered budget so that mass stays small relative to the estimated hidden-set size.
- **Pre-registered variant selection.** Three field variants are computed once: `euler_only` (H4-A0), `euler_x_terrain` (H4-A1, terrain concurrence floor 0.35) and `euler_x_terrain_x_gap` (H4-A2, with mapping-gap weighting). Selection rule fixed now: rank by pooled catalogue-block DTI; if the top two are within 0.005, prefer the geologically-motivated concurrence variant (H4-A1 or H4-A2) over `euler_only`. This is a proxy-instrument selection and is reported as such.
- **Expected improvement / cost.** Directionally positive versus a pure single-field Euler cloud because three independent evidence classes must agree before a dot is emitted; magnitude unknown. Cost: medium (Euler solves ~15 s per configuration on this CPU; lineament transform and audits are minutes).

### H4-B (rank 2): mapping-completeness (catalogue-gap) weighting of the contact field

- **Layers.** Official USGS QFFD scale split (`coarse_trace`, `fine_trace`, `lower_certainty_trace` rasters) crossed with the H4-A field.
- **Physical signature / rationale.** Fault-database completeness is spatially non-uniform: the footprint carries 58,251 QFFD coarse-scale px but only 2,754 fine-scale px, i.e. this region is compiled at reconnaissance scale. Areas of dense coarse mapping are "already found"; areas where strong structural evidence exists but the compile is silent are the gaps the hidden set is drawn from.
- **Difference from this repo.** The repository used an SGMC *gap* as a truth proxy (circular). H4-B uses the QFFD scale split as a *covariate on candidate emission*, never as truth.
- **Cost.** Low (rasters already available); risk: QFFD density may track real fault density rather than mapping effort.

### H4-C (rank 3): multi-physics alteration concurrence (Euler ⊕ radiometric ratios ⊕ upward-continued TMI)

- **Layers.** GeoDAWN airborne radiometrics (K, Th, U, TC) and contractor ratio grids (Th/K, U/K, U/Th, TMI upward-continued 150 m) — official USGS data release DOI 10.5066/P93LGLVQ.
- **Signature.** Along-fault potassium enrichment / uranium depletion relative to thorium (clay and silica alteration) forms linear anomalies; requiring an Euler contact to coincide with a ratio anomaly targets structures with hydrothermal alteration, which is what makes a fault geothermally relevant.
- **Difference.** The sibling repository 15GEMSDOE ran a `conj_alteration_mag` conjunction (0.0782 live); this repository has never implemented an alteration term, and the H4-C variant registers the ratio *lineament* rather than the raw ratio value.
- **Cost.** Medium; the radiometric grids are already local.

### H4-D (rank 4 — blocked, not viable in this sandbox): focal-mechanism / seismicity lineament transform

- **Layers.** USGS ComCat origins and focal mechanisms plus challenge bands `deq_n100a15`, `ieq_n100a15`.
- **Signature.** Orientations of nodal planes, accumulated as a directed lineament vote, sharpen active-fault geometry beyond the smooth earthquake-density bands.
- **Viability check performed 2026-10-06:** `earthquake.usgs.gov` returns HTTP 000 from this sandbox (egress-restricted), and only the count endpoints were reachable through the research web tool in the earlier session; the H3 coverage gate (≥50 usable mechanisms spread over ≥8 of 24 blocks) is still unverified. The specific free official source required is the USGS ComCat FDSN event web service `https://earthquake.usgs.gov/fdsnws/event/1/`. **Not viable here without an unrestricted runner.**

### H4-E (rank 5): blind basin-interior fault detection

- **Layers.** `depth_to_base_surf` (15), `cond_surf` (17), `iso_grav_anom_hg` (18), `det_elev` (12).
- **Signature.** A step/flexure in the basement surface beneath basin fill, aligned with a conductivity discontinuity and a gravity gradient, is the classic signature of a buried normal fault — the class of structure a surface compilation is least likely to contain.
- **Difference.** No repository arm has ever used the sediment-thickness or conductivity bands.
- **Cost.** Medium–high; these layers are smooth and the `tc`/`cond_surf` semantics in the band inventory are partly inferred, so the arm needs its own semantics audit first.

**Ranking by expected DTI improvement per unit cost:** H4-A > H4-B > H4-C > H4-E > H4-D (blocked). Only H4-A is implemented in this session; H4-B is applied as a weighting variant inside H4-A (pre-registered above); H4-C is left as a registered, unimplemented candidate with its data already verified present.

### Append-only amendment — H4-A result, attribution integrity, and the promotion-gate blocker (2026-10-06 UTC)

**1. Binding integrity rule.** A leaderboard row establishes an *account-level* public score only. It does not identify a raster unless the organizers link the submission file. Any filename↔score pairing without such a link is an **unauthenticated attribution** and may not be used to choose emission rules, weights, thresholds or field combinations. `docs/data/feed-20261005.json` already records this and specifically forbids attributing the 0.2778 row to `gemsdoe32-h33-h33-2-b2-…`.

**2. Retraction.** This session ran an emission calibration over ten such unauthenticated pairs (`scripts/audit_unauthenticated_attributions.py`, output `docs/data/unauthenticated_attribution_audit.json`) which produced strong apparent correlations (median distance to catalogue rho +0.988; mean detrended elevation +0.988; emitted pixel count −0.976; fraction of dots within 300 m of the catalogue −0.988). **Withdrawn as evidence.** Both sides of those correlations are unauthenticated, so they establish neither a design rule for emission nor any statement about the registered instruments (in particular, they do not show that the catalogue-block instrument is anti-predictive, because the same unauthenticated labels appear on the other side of that comparison). The weights module derived from them (`src/gemsdoe40/h4_weights.py`) was deleted, and the design (a "match the known-best files" emitter, drafted as H4-F) was abandoned **before any raster was written**. The analysis is retained only as a documented record of a prohibited method.

**3. H4-A result — negative, closed.** Pre-registered comparator: 4×6 contiguous blocks with a three-cell guard, pooled official DTI on `data/labels.tif`. At the ladder-maximum budget of 90,000 dots: A0 `euler_only` **0.0172**, A1 `euler_x_terrain` **0.0169**, A2 `euler_x_terrain_x_gap` **0.0168**; the registered tie rule fired (top two within 0.005) and A1 was emitted. Same-mass random control **0.1257**, gradient-top-K baseline 0.0365, contiguous (non-NMS) top-90k extraction 0.0372; 2-fold held-out block selection 0.0106–0.0107. The detector is **7.4× worse than chance** on its own comparator, so H4-A is closed and **not promotable**; its artifact was moved out of `docs/downloads/` to `work/withdrawn/` so it cannot be mistaken for a shippable candidate. Diagnosis: the KDE field is positive on 2,255,615 px (43.6 % of the footprint) and separates catalogue-adjacent cells from the rest by under 4 px of mean distance (AUC 0.534 against the catalogue), so the field is not spatially selective; the failure is upstream of the emitter.

**4. Promotion gate is currently un-runnable.** The frozen gate's proxy raster (owner-derived SGMC mirror, SHA-256 `26d142c4c93282cd94f6950ab96f22aeff59fbbea523d43d662e76fa1b161b5c`) lived in `/tmp/proxy-data/`, which does not persist between sandbox turns, and the source host is outside the sandbox egress allowlist. The surviving SGMC raster (`work/ext/derived_sgmc_faults_100m_u8.tif`, SHA-256 prefix `643cbe99…`) is the circular one used to build the frozen incumbent, so scoring against it cannot promote anything. Until the pinned proxy is restored or an organizer-linked instrument exists, **no arm may be promoted to a submission on local evidence**; diagnostic rankings computed on the surviving raster are recorded in `docs/data/proxy_ranking.json` with that caveat attached.

**5. H4 slate status.** H4-A closed (negative). H4-B retained only as a covariate (its gap weighting is inside H4-A and did not help). H4-C (radiometric alteration concurrency) still registered and unimplemented, data verified present. H4-D blocked (egress). H4-E registered and unimplemented. H4-F abandoned unbuilt, and its motivating correlations withdrawn under rule 1.

**6. Organizer rulings recorded this session (context for Phase 2, no new arm).** Two DrivenData staff statements, both from the official forum, change how a Phase 1 submission should be composed and are recorded here so later arms do not re-derive them:

- *Definition of a new fault.* Staff (chrisk-dd, topic 11536, 2026-09-23): "new fault" means "any fault pixel not already captured by USGS/INGENIOUS" and **can include newly mapped geometry of an existing fault system** — i.e. continuations, splays and parallel strands of mapped faults are admissible targets, they are simply not scored where they coincide with the masked known-fault pixels.
- *Phase 2 uses expert review of Phase 1 submissions.* Staff (chrisk-dd, topic 11527 post 7, 2026-09-23): "We're not sharing details about the data sources, fault types, or coverage behind the test faults beyond what's in the problem description", and "the largest prize pool (Phase 2) will use a test set that is updated by expert review of all Phase 1 submissions, so your fault predictions have an impact on final evaluation even if they are not the most performant in Phase 1."

Consequence for arm design: a Phase 1 submission is also a *proposal set* for expert review. Mass that is spatially diffuse, or that cannot be defended from surface or geophysical evidence, adds Phase 2 review burden without adding admissible targets, while coherent, evidence-backed candidate traces (even at modest Phase 1 score) can be adopted into the Phase 2 label set. This is a strategic judgement, not a measured effect, and it does not relax any registered gate.

**7. Correction to item 4 (corpus ranking not completed).** The equal-mass ranking of all 279 cached priors on the surviving SGMC raster was started (`scripts/rank_priors_proxy_instrument.py`, 45,000-dot equal-mass emission, 4x6 blocks with guard 3) but did not complete inside the sandbox CPU budget: 60 of 279 rasters were scored before the run was stopped. The partial record is `work/rank_probe.json` plus the process log (no committed artifact), and **`docs/data/proxy_ranking.json` therefore does not exist** — item 4's reference to it is superseded by this correction. The partial result already shows what item 4 predicted: the top of the list is the incumbent's own lineage (`GEMSDOE3 gapfinder-v2-fusion`, DTI 0.4174 against a three-draw same-mass random control of 0.0825, i.e. 5.1x), which is circular by construction. Re-run the script on an unrestricted runner before quoting any corpus-wide ranking.


## 2026-10-06 — contact-offset H4–H7 slate (separate from H4-A and H7 gravity-context)


A separate frozen [H4–H7 registration](h4-preregistration-20261006.md) preceded the new implementation. It ranks offset-aware rank-adaptive contact Euler, finite gravity-step inversion, depth-cloud plane geometry, and TMI/RTP representation stability. Only H4 was implemented. H6's ground-surface projection remains blocked pending verified survey datum.

H4 now exists as a genuinely distinct, format-valid research TIFF, with byte-identical independent regeneration. It **fails** blocked proxy promotion; it is not a candidate for a weekly upload. See the [negative result](../reports/h4-results-20261006.md). No post-score hyperparameter search or H4 re-emission was performed.

### Integration correction — 2026-10-06

The **contact-offset H4** in `h4-preregistration-20261006.md` is TMI plus dG/dz with free A and rank adaptation. It is distinct from the concurrent **H4-A** terrain/QFFD contact-network arm and the merged **H7 gravity-context** arm. H7 in the contact-offset slate means unimplemented TMI/RTP representation stability, not the already implemented gravity-context experiment. Preserve all records; do not mix their scores or labels.

The formerly unavailable `26d142c4…` SGMC proxy has now been autonomously recovered from the pinned public GEMSDOE30 commit; see `../data/acquisition-20261006.json`. Its frozen 0.8359066541 incumbent score is reproduced. This resolves missing data, not circularity: the incumbent uses SGMC, and only 16 blocks contain truth while the rule requires 18 wins. The instrument is retained for historical diagnostics only and is **retired for promotion**. No replacement gate based on unauthenticated owner file/score pairs is adopted.

## 2026-10-06 session 2 — H8–H13 slate

The session-2 slate (H8 trace-locked Euler depth consensus, H9 blind basement-flexure, H10 geodetic strain lineaments, H11 seismic corridor × Euler, H12 ComCat mechanisms) plus the H13 value-concentrated crest emission and the CAT-HID instrument repair are registered in [`h8-preregistration-20261006.md`](h8-preregistration-20261006.md). Results: both implemented candidates HOLD at the frozen proxy gate; both pass the repaired CAT-HID skill check. See [`../reports/session2-results-20261006.md`](../reports/session2-results-20261006.md).
