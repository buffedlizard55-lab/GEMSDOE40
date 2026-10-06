# Preregistered hypotheses — 6 October 2026 (second sitting)

Frozen before implementation. Written for the standing brief that requires 3–5 candidate
geological hypotheses, each naming its layers, the physical signature it targets, why it should
catch a fault **missing** from the USGS/INGENIOUS catalogue rather than one already in it, and
how it differs from everything already tried in this repository — then ranked by expected
improvement per unit cost, with the top candidate validated on the spatially blocked holdout
**before** any weekly submission slot is spent.

**Validation rule used.** A candidate must beat the current holdout comparator on the frozen
4 × 6 blocked instrument with a 3-pixel guard, *and* the instrument used to project a live score
must reach leave-one-out Spearman ≥ 0.80, before a slot is contemplated. The second condition
is not met by any instrument in this repository (best LOO 0.705), so every candidate carries
**HOLD** until the owner decides otherwise.

**Ranking is judgement; the measured column is separate.** Nothing below is a leaderboard
forecast.

---

## H41 — IMPLEMENTED (this release) · cost: low (hours)

- **Layers.** Euler cloud from `rtp` (band 2), `tmi` (14), `iso_grav_anom` (13); ranker features
  from `tmi_hg` (3), `geod_2ndinv` (4), `iso_grav_anom_slope` (5), `tc` (6), `iso_grav_anom_vg`
  (11), `det_elev` (12), `depth_to_base_surf` (15), two cross-field gradient ratios, Euler
  multi-window summaries and DEM gradient/Laplace terms (27 total).
- **Signature.** Depth-labelled contact source positions from **SI = 0 Euler deconvolution**
  (Reid, Allsop, Granser, Millett & Somerton, *Geophysics* 55(1), 1990), clustered by shallow,
  mutually consistent, multi-scale-stable kernel density; emitted as a sparse binary dot set.
- **Why it catches an unmapped fault.** A fault under cover has no surface scarp but still
  produces a coherent potential-field contact line; the provided catalogue is surface-mapped, so
  such a contact can fall entirely outside it. Requiring a solution to survive three window
  scales suppresses single-window artefacts that dominate naive Euler fields.
- **How it differs from everything already here.** Earlier arms either threshold a surface
  derivative or cluster Euler solutions without a depth-consistency requirement; H41 requires
  three-scale stability and ranks the density against public SGMC faults that are missing from
  the provided catalogue, out of fold (pooled OOF AUC 0.7488).
- **Measured result (receipts, not promises).** Blocked proxy mean **0.1366** against 0.0902 for
  the owner's best scored file (H33-B2), 0.1070 for H27-4 and 0.1395 for H40-E; 8/16 strict
  block wins; proxy DTI 0.1284 at 40,000 emitted dots; unique versus all 51 staged priors
  (max |r| 0.0164, max support Jaccard 0.0119).
- **Frozen before the run.** Windows 10/16/24 px, stride 4 px, SI = 0, depth ≤ 1500 m, 3-pixel
  emission spacing, catalogue flank 1.0 px. The flank decision (1 px versus 2 px) was taken
  *before* the publish run; both were measured (H41 0.1272 on 5/16 blocks, H41c 0.1298 on 7/16).
- **Receipts.** `docs/data/h41-generation.json`, `docs/data/h41-audit.json`,
  `docs/research/h41-preregistration-20261006.md`, `docs/reports/h41-results-20261006.md`.

## H42 — REGISTERED · blocked here · cost: high (new survey metadata)

- **Layers.** Euler clouds + the official airborne-magnetic/radiometric survey drape geometry.
- **Signature.** Dip azimuth and dip from the horizontal drift of Euler solutions with solved
  depth, projected to a surface trace.
- **Why it catches an unmapped fault.** A dipping blind fault migrates laterally with depth;
  projecting the fitted contact plane to the surface places dots where no scarp exists yet —
  exactly the population an expert unmapped-fault set contains.
- **How it differs.** No arm in this repository fits orientation: the rank-aware solver discards
  the unidentifiable along-strike coordinate instead of projecting it.
- **Blocked by.** The official GeoDAWN geometry/per-flight-line metadata
  (<https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and>).
  Obtainability was checked in this environment: every USGS host returns no response, so the
  idea is **not** called viable here. The source itself is free and official; a session with
  USGS access should re-test obtainability first.

## H43 — REGISTERED (cheapest next test) · cost: low (bands already in the grid)

- **Layers.** `tmi_vg` (9) and `iso_grav_anom_vg` (11) — vertical gradients of the magnetic and
  gravity fields.
- **Signature.** Rock-property contrast: the ratio of magnetic to gravity vertical gradient
  separates hydrothermal alteration zones from unaltered intrusive contacts.
- **Why it catches an unmapped fault.** It removes the dominant false-positive class of shallow
  contact solutions (intrusive contacts), so surviving dots concentrate on fault-like
  boundaries rather than on every lithologic edge.
- **How it differs.** H41 uses both bands as independent features; the ratio is a physical
  discriminant rather than another correlated channel.

## H44 — REGISTERED · cost: low (needs a DEM list) · partially blocked

- **Layers.** `labels.tif` + `det_elev` (12) + Euler clouds.
- **Signature.** Fault-network topology — tips, step-overs and along-strike gaps — with a DEM
  lineament required to bridge the gap.
- **Why it catches an unmapped fault.** Unmapped faults are frequently along-strike extensions of
  mapped ones; the owner's own measurement shows dots inside 200 m of the catalogue hurt, so the
  extension must be geometrically constrained, not assumed.
- **How it differs.** The repository's earlier gap-closure arm used topography alone (0.2449
  live); H44 requires an Euler depth cluster to bridge the gap before a dot is emitted.
- **Blocked by.** The organizer's 1 m DEM link list lives behind the login-walled data page; a
  public equivalent (e.g. 3DEP) can substitute but was not verified from this sandbox.

## H45 — REGISTERED · cost: low (bands already present)

- **Layers.** `deq_n100a15` (10) and `ieq_n100a15` (16) — earthquake-density bands.
- **Signature.** Seismicity lineaments tested against a Poisson null on the same footprint.
- **Why it catches an unmapped fault.** Active structures can post-date the Quaternary
  compilation used for the labels; the expert set may include faults with instrumental
  seismicity and no mapped scarp.
- **How it differs.** No arm in this repository has ever ranked on the seismic-density layers.

---

## Decision recorded

H41 was implemented, audited and published as the current downloadable release; the blocker is
that no local instrument clears the 0.80 projection bar, and a prior unscored Euler arm (H40-E,
0.1395 blocked mean) still leads H41 on this repository's own proxy. Per **Maximize P(Win)** and
**Own the Outcome**, no weekly slot was spent and no organizer score is claimed. H43 is the
cheapest next test; H42 is the highest-upside idea and is blocked by a named, free, official
source that was verified unreachable here.
