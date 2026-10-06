# Frozen H8 lineament implementation and the ranked candidates — 2026-10-06 UTC

> Session-2 sibling register: [`h8-preregistration-20261006.md`](h8-preregistration-20261006.md) (retained,
> byte-identical, with its own generation-time hashes). This document registers the **lineament-weighted
> cross-family SI = 0 depth-clustering** variant implemented in this session, its frozen settings and its
> promotion rules, before the run.

**Numbering.** The session-2 register already used H8–H13 (H9–H11 there were closed NEGATIVE at frozen
pre-screens, H13 is the crest-binary artifact). The four proposals below are therefore numbered **H14–H17**
to avoid any collision; the artifact this document registers keeps its descriptive filename
`…h8-euler-lineament-depthcluster…` because it is a variant of the same H8 Euler depth-clustering family.

**Registration:** 2026-10-06 (UTC). **Rule applied:** *Maximize P(Win)*, *Own the
Outcome*. H8's generator settings were frozen in
`src/gemsdoe40/h8_euler.py` / `scripts/run_h8_lineament.py` **before** the run;
the run receipt `work/h8lineament/h8-lineament-receipt.json` records them
verbatim and the candidate is written only into ignored `work/`. H8 reads no
prior prediction, no holdout raster and no leaderboard value. Nothing in this
document is a competition score.

## 0. What was measured before any hypothesis was written (this session)

| measurement | result | where |
|---|---|---|
| Frozen H4 SI=0 contact cloud vs the mapped catalogue | 46,656 solutions; weighted mean distance **32.95 px**; P(d < 100 m) **0.0099** vs random-footprint **0.0115** — at chance | `work/diag_h4_cloud.py` |
| Same solver, 600 × 550 cell crop around the largest mapped fault, RTP / TMI / isostatic gravity, windows 9/15/25 | P(<100 m) 0.2–1.3 % against a random baseline of 0.8–1.5 %; P(<300 m) 2.0–8.7 % against 6.6–8.3 % — **no window size or field beats chance** | `work/euler_probe.py` |
| Control test on an exact analytic SI=0 contact field, `T = A·atan2(x−x₀, z−z₀)` | depth recovered exactly (500 m at every window); position recovered to the noise floor; the solver is not broken | `work/euler_control.py` |
| Transform skill against the catalogue, 19 bands × 6 transforms | best per-dot credit is `geod_2ndinv raw_abs` (w = 0.168 at 20 k, hit-rate 5.1 % vs 1.15 % random), a **regional** concentration effect; every lineament transform (gradient magnitude, Laplacian, tilt) sits at w ≈ 0.04–0.06 | `work/transform_skill.py` → `work/transform_skill.json` |
| Ranking power of four instrument families on the 16 recorded scores | catalogue credit ρ = −0.285; SGMC-off-catalogue credit ρ = +0.491; mass ρ = −0.650; score-anchored peer consensus (LOO) ρ = +0.54…+0.70 — **none reaches the ρ ≥ 0.8 the family's own audit requires** | `work/prior_catalogue_matched.json`, `work/truth_model.json` |
| Implied hidden credit per emitted pixel (metric inverted on the recorded scores, N = 5,800) | best recorded file **0.090**, worst **0.016**; beating 0.3195 at the same mass needs **0.103** | `docs/research/h8-analysis-20261006.md` §4–5 |

## 1. Ranked hypotheses

Ranking is by expected DTI improvement per unit of implementation cost. The
"layer" names are the exact `band_name` tags of `data/training_features.tif`
(19 bands, verified this session); no band number is assumed anywhere.

### Rank 1 — H8 (implemented today): lineament-weighted, cross-family-corroborated SI=0 Euler depth-clustering

* **Layers:** `rtp` (reduced-to-pole total magnetic intensity) and
  `iso_grav_anom` (isostatic gravity anomaly) differentiated once in the Fourier
  domain by the `+|k|` operator, i.e. the contact-like gravity arm; `tmi` and
  the supplied `iso_grav_anom_vg` remain independent cross-checks.
* **Physical signature:** depth-labelled *contact* Euler solutions
  (SI = 0; Reid, Allsop, Granser, Millett & Somerton 1990, eq. 2 including the
  arbitrary contact offset A) — the structural index the 1990 paper uses for a
  fault-like contact. Not an edge map and not gradient thresholding: the
  output is a solution cloud, then a continuous kernel-density field.
* **Why a *missing* fault and not a catalogue fault:** a potential-field contact
  is a subsurface boundary; the catalogue was compiled from surface mapping. A
  fault under basin fill has no surface trace but still bounds a magnetisation
  or density contrast. The emission additionally excludes the pixel-exact
  catalogue mask and its immediate neighbours, because the metric cannot award
  credit there (DrivenData staff, 2026-09-16).
* **Difference from everything already in this repository:** H4 and H7 run the
  same solver but emit the full continuous KDE and treat all solutions alike.
  H8 is the first artifact that (i) weights each solution by the *geometry* of
  the cloud (`1 − λ₂/λ₁` lineament coherence of its 500 m neighbourhood),
  (ii) halves the weight of solutions the other physics family does not
  independently corroborate within 400 m and 600 m depth, (iii) smears the
  weighted cloud with an *anisotropic* kernel aligned to the local lineament,
  so solution trains become trace-like ribbons, and (iv) emits the metric-tuned
  sparse support (2.8 px maximum-separation dots, fixed 40,000 budget) with
  values that are the normalised density (continuous, [0, 1]).
* **Expected improvement: low–moderate, stated honestly.** The probe in §0 says
  this solver's solutions land on mapped faults at chance in this region at
  every window tested. H8 is the mandated construction and a genuinely new
  pattern; it is not a forecast.
* **Cost:** paid (medium).

### Rank 2 — H14: tip, step-over and along-strike *extension* of mapped systems

* **Layers:** the provided catalogue geometry (`data/labels.tif`) as the seed,
  with `det_elev_slope`, `rtp` and the gravity-gradient ridges as the geometric
  corroboration along the projected strike.
* **Physical signature:** pure geometry — strike-projected tip extensions of
  1–3 km, relay/step-over corridors between overlapping strands, and
  along-strike continuation of short mapped segments.
* **Why this is exactly what the hidden set contains:** DrivenData staff,
  2026-09-23, thread 11536 ("Where do you draw the line?"), post 2: "For the
  purposes of this competition, 'new fault' means 'any fault pixel not already
  captured by USGS/INGENIOUS' and can include newly mapped geometry of an
  existing fault system." (<https://community.drivendata.org/t/where-do-you-draw-the-line/11536/2>) The catalogue pixels themselves are masked
  out of scoring, so an extension beyond a mapped tip is scored normally while
  the mapped trace costs nothing.
* **Difference from prior work:** the family's best file (H33-B2) *deleted* all
  mass within 200 m of the catalogue. No artifact in the corpus emits
  tip-extension or step-over geometry. This is the deliberate opposite sign of
  the same measurement.
* **Validation available now, with real truth:** leave-the-tips-out — cut the
  terminal 25 % of every catalogue trace in a spatially blocked fashion, build
  the extension field from the truncated geometry only, and score the withheld
  tips with the official metric. No new data needed.
* **Cost:** moderate. **Expected improvement: highest of the five.**

### Rank 3 — H15: finite-step gravity inversion at basin margins

* **Layers:** `iso_grav_anom`, `iso_grav_anom_vg`, `iso_grav_anom_hg`,
  `depth_to_base_surf`.
* **Physical signature:** the *finite* density step of a fault block, not the
  infinite contact. Reid & Thurston (2014) show the finite step has SI = −1 and
  needs the generalised multi-edge treatment; the SI = 0 vertical-gradient arm
  used by H4/H7/H8 is only a local top-edge approximation, which is a known
  weakness of the current line.
* **Why missing from the catalogue:** basin-margin faults under fill are the
  class that surface compilation systematically under-maps; the
  `depth_to_base_surf` band gives an independent geometric prior on where the
  step should be.
* **Difference:** every Euler artifact in the repository assumes the
  infinite-contact approximation. H15 is a different inverse problem.
* **Cost:** high. **Expected improvement: moderate.**

### Rank 4 — H16: 1 m lidar scarp re-mapping over the full footprint

* **Layers:** 3DEP 1 m DEM (competition `1m_DEM_links.csv`, login-walled) or the
  public USGS 3DEP service; the repository's own `lidar_scarp` feature family
  (`upface_max` AUC 0.5914) is the prior partial implementation, covering
  26–32 % of the footprint.
* **Physical signature:** up-face/down-face scarp asymmetry, cross-scarp
  curvature and differential relief at 1 m, aggregated to the 100 m grid — the
  direct surface expression of the fault type the experts label.
* **Why missing from the catalogue:** lidar-visible scarps in young deposits
  are precisely the recent faults that a pre-lidar compilation missed.
* **Difference:** the challenge bands `det_elev` / `det_elev_slope` are coarse
  regional rasters (catalogue AUC 0.4855 / 0.5670); no artifact in the corpus
  performs a scarp re-mapping over the whole footprint.
* **External data source, named and checked:** USGS 3DEP 1 m DEM, free and
  official — <https://www.usgs.gov/3d-elevation-program> and the TNM Access API
  <https://tnmaccess.nationalmap.gov/api/v1/products>. **Obtainability from
  this sandbox: not possible** — direct HTTPS to `prd-tnm.s3.amazonaws.com`
  fails at the TLS handshake here (recorded in the session log); the data is
  obtainable on an unrestricted machine, and the competition's own data tab
  ships the CSV of pre-cut DEM links. Do not claim it was downloaded here.
* **Cost:** high (large data volume). **Expected improvement: high if the
  hidden labels are scarp-derived, zero otherwise.**

### Rank 5 — H17: microseismicity-aligned structures

* **Layers:** `ieq_n100a15` (smoothed earthquake intensity) and `deq_n100a15`
  (distance to the nearest event), with `geod_shearrate` as a secondary field.
* **Physical signature:** a fault that is active but unmapped concentrates
  microseismicity; the intensity ridge, not the raw event count, is the signal.
* **Why missing from the catalogue:** the USGS Quaternary compilation is
  geomorphology-led; an active structure under young fill can be absent from it
  while still generating events.
* **Difference:** no artifact in the corpus uses the seismicity bands as a
  *primary* discriminator (measured catalogue skill is weak: `ieq_n100a15`
  gradient magnitude w = 0.0531 at 20 k).
* **Cost:** low. **Expected improvement: low.**

## 2. Frozen H8 settings (recorded before the run)

`src/gemsdoe40/h8_euler.py`: families `rtp` (continuation 150 m, order 0,
windows 9/15/21, max depth 3,000 m) and `iso_grav_anom` (continuation 500 m,
order 1, windows 15/21/31, max depth 5,000 m); SI = 0 with contact offset;
stride 4; FFT pad 192; depth decay 1,200 m; cluster radius 300 m with a depth
tolerance of max(200 m, 0.35·depth) and ≥ 3 neighbours including one from a
different window; coherence radius 500 m with ≥ 5 neighbours; corroboration
radius 400 m XY and 600 m depth; uncorroborated solutions × 0.5; anisotropic
kernel σ = 1.6·(1 + 4·coherence) px along the local lineament axis and 0.55 px
across it; maximum-separation emission spacing 2.8 px; fixed budget 40,000
dots; values = density / median(density of accepted dots), clipped to [0, 1];
exclusion = the pixel-exact catalogue mask dilated by 1 px.

## 3. Promotion rules applied (unchanged from the repository's standing rules)

1. Generation writes only into ignored `work/`; no prior prediction, holdout or
   score enters construction (`construction_reads_labels_or_priors: false`).
2. Publication into `docs/downloads/` requires the format gate (single float32
   band, EPSG:32611, 3730 × 3292, exact sample transform, finite [0, 1] at all
   5,167,373 footprint cells, NaN outside) **and** the novelty gate (no exact
   duplicate; max |Pearson| ≤ 0.95; max top-k Jaccard ≤ 0.50) against every
   cached prior raster.
3. No local instrument currently available to this programme reaches the
   ρ ≥ 0.8 ranking power the family's audit demands, so no local number is
   described as evidence of a leaderboard improvement. See
   `docs/research/h8-analysis-20261006.md` §6.
4. The submission slot is a scarce resource: the generated file is published
   for inspection and download with an explicit status, not as a forecast.
