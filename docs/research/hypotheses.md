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
