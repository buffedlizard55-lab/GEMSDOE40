# GEMSDOE40 — Euler depth-cluster fault prediction

> **Read this README and both retained briefs at the beginning of every session.** See [AGENTS.md](AGENTS.md). Core Values: **Maximize P(Win)** and **Own the Outcome**.

## Download the submission GeoTIFF

**[Download H41 — the unique, format-valid submission GeoTIFF](docs/downloads/h41-msst-gated27-923c57ab-zeros.tif)** · [Live project site](https://buffedlizard55-lab.github.io/GEMSDOE40/docs/index.html) · [Executive summary / how to submit](docs/executive-summary.html) · [NaN-outside twin](docs/downloads/h41-msst-gated27-923c57ab-nan.tif) · [.zip containing one GeoTIFF](docs/downloads/h41-msst-gated27-923c57ab-zeros.zip)

**Unverified, but the strongest construction this arm has produced — and the repository's HOLD gate still applies.** The file is unique against every retrievable prior output, passes every format check on its published bytes, and improves on the owner's best scored file (0.2778 → see §3) on the two instruments this repository can recompute. It does **not** beat one prior, never-scored Euler arm of a sibling repository (H40-E) on the blocked proxy instrument, and no locally reproducible instrument reaches the pre-registered 0.80 leave-one-out rank correlation that this repository requires before it certifies a slot. No weekly slot was used here; no organizer score is known for any file named in this repository.

| Item | Measured result |
|---|---|
| File | `docs/downloads/h41-msst-gated27-923c57ab-zeros.tif` (204,000 bytes) |
| Submission name | `GEMSDOE40-H41-MSST-EULER-923c57ab` |
| Short note (163 chars) | `GEMSDOE40 H41 | SI0 contact Euler depth-clusters (3 windows, 356,650 solutions); blocked-OOF ranker; 40,000 dots, 3 px spacing; unique vs 51 priors (max |r| 0.014)` |
| Format | Single-band float32, EPSG:32611, 100 m, 3,730 rows × 3,292 columns, transform `(100, 0, 243350, 0, -100, 4508550)`; **every** cell finite in [0, 1] (40,000 cells = 1.0, everything else 0.0 inside the footprint, 0.0 outside) |
| File SHA-256 | `1a64839d5195728df53303db5d6e9132efa2aef863c2416a41e53cc6eadd74d5` |
| Euler deconvolution | SI = 0 contact form (offset A, Reid et al. 1990 eq. 2 / step 3b) on bands 2 (rtp), 14 (tmi), 13 (iso_grav_anom); windows 10/16/24 px, stride 4 px, gates: analytic-signal percentile 72, relative depth error ≤ 0.22, depth 80–2,200 m; **356,650** solutions |
| Depth clustering | Per-window Gaussian KDE (σ = 1.7 px) with shallowness × quality × tightness × depth-consistency weights and magnetic/gravity concordance; combined by **geometric mean across all three windows** (2,187,859 cells supported at every scale); solutions ≤ 1,500 m only |
| Ranking | Spatially blocked (5 × 6) out-of-fold L2 logistic discriminant, 27 features, trained on public USGS SGMC faults **absent from the provided catalogue**; pooled out-of-fold **AUC 0.7488** (mean of blocks 0.7615) |
| Emission | Value-ranked Poisson-disk thinning at 3 px spacing (the metric's own R = 300 m) inside the Euler-gated support (2,101,752 cells, catalogue flank 1 px); mass chosen by the recalibrated instrument: **40,000 dots**, w = 0.0987 |
| Blocked proxy comparison (4 × 6, 3 px guard) | H41 **0.1366** mean · H40-E 0.1395 · H27-4 0.1070 · H33-B2 0.0902 (7 of 16 truth-bearing blocks strictly beat all three comparators) |
| Whole-map proxy DTI | H41 **0.1284** · H40-E 0.1466 · H27-4 0.1069 · H33-B2 0.0913 |
| Instrument estimate | H41 **0.2282** · H33-B2 (owner-reported 0.2778) 0.1820 · H27-4 (0.2708) 0.1999 — recalibrated saturating model, leave-one-out Spearman **0.705**, RMSE 0.054 |
| Novelty | Maximum |Pearson| **0.0164**, maximum support Jaccard 0.0119, maximum top-budget Jaccard 0.0097 against 51 same-grid prior rasters; `is_new = True` |
| Decision | **HOLD on this repository's gate** (instrument below the 0.80 bar; H40-E leads the blocked proxy). Publish → the owner decides. `slot_eligible: false`, `organizer_score: null` |

### 1. Why the previous best (0.2778) was the best — measured, not guessed

The published metric reduces to `DTI = T / (0.2T + 0.2(M − C) + 0.8N)` with `T` the distance-weighted true-positive credit, `M` the emitted mass, `C` the credit inside that mass and `N` the hidden truth size. Four measured facts explain the ordering of the 16 owner-reported scores in this repository's anchor corpus ([receipt](docs/data/live-transfer.json)):

1. **Mass dominates.** Spearman correlation between emitted mass and live score is strongly negative; the measured live/proxy-DTI ratio falls monotonically from 2.53 (40 k dots) to 0.28 (146 k dots). The best file emits only 37,654 pixels.
2. **Arrangement is worth ~0.2.** At the *same* 37,654-pixel budget, a scattered control scores 0.0778 while the structure-aligned file scores 0.2778 ([measured](docs/data/h33-measured-analysis.json)).
3. **The best file adds nothing.** It is the 40,199-pixel base with 2,545 dots within 200 m of the known catalogue deleted — the minimum remaining distance to a mapped fault is 223.6 m.
4. **The public surrogate is not a promotion instrument.** The whole-map SGMC-off-catalogue DTI ranks the scored anchors *backwards* out of fold (LOO Spearman −0.897); at fixed mass the per-dot credit `w` does rank them (Spearman +0.705). Optimising the surrogate therefore cannot justify a slot, which is why this repository publishes the number and the error bar instead of a promise.

**Can 0.2778 be beaten?** Only by more hidden-truth credit per unit mass. H41 does exactly that on paper — same mass as the best file with a 21 % higher per-dot credit `w` (0.0987 vs 0.0816 at 40,199 dots), hence an instrument reading of 0.2282 versus 0.1999 for the same-family 0.2708 file and 0.1820 for the 0.2778 file. The instrument under-predicts the two best anchors by ~0.09, so if that bias is systematic H41's live score could sit in the low 0.3s; if it is not, H41 lands near the owner's best. That is the honest uncertainty: **one extrapolation, clearly labelled, not a forecast.**

### 2. What is verified and what is not

* Verified on the published bytes: format legality ([audit](docs/data/h41-audit.json)), uniqueness against all 51 retrievable prior rasters, the blocked proxy comparison, the metric algebra, and the instrument's own leave-one-out performance.
* Not verified: any organizer score. No competition account, no submission and no scraping exists in this environment, and DrivenData's [terms of use](https://www.drivendata.org/termsofuse/) forbid automated monitoring.
* Known deficiency of the passed gate: the pre-registered promotion bar (instrument LOO Spearman ≥ 0.80) is **not** met at 0.705, and H40-E — a prior, never-scored Euler arm from the sibling GEMSDOE39 repository — still leads the blocked proxy mean (0.1395 vs 0.1366). Both are published rather than smoothed over.

### 3. The submission-format bug, answered

The portal error "Predicted values must be in range [0, 1]" has two verified causes in this repository: (a) `training_features.tif` stores the sentinel **−3.4028235e+38** as data outside its valid area, so any raster assembled from it wholesale contains values near −10³⁸; and (b) **NaN fails a numeric range test** (`0 <= nan <= 1` is false), so a validator that tests the whole array rejects the sample submission's own NaN-outside encoding. H41 writes every cell explicitly: exactly 1.0 for the 40,000 emitted dots and exactly 0.0 everywhere else inside the footprint — no NaN, no sentinel, no value outside [0, 1]. A NaN-outside twin is provided for completeness.

### 4. Ranked hypotheses (see [hypotheses page](docs/hypotheses.html) and [pre-registration](docs/research/h41-preregistration-20261006.md))

**H41** (implemented, this release) multi-scale-stable SI = 0 depth-cluster emission → **H43** cross-field vertical-gradient ratio as a physical discriminant (low cost, bands already present) → **H45** seismic-density lineaments against a Poisson null (low cost) → **H44** catalogue gap closure *bridged by* an Euler depth cluster (low cost, needs the pinned 1 m DEM list) → **H42** dip projection from the depth-migration of Euler solutions (high cost; blocked here by GeoDAWN survey-geometry metadata that is named and linked but unreachable from this sandbox).

## Concurrent sibling release — H45 (merged 2026-10-06, preserved verbatim)

The H45 arm landed on `main` while this release was in review. It is a separate session's file with
its own name, bytes, receipts and HOLD decision; it is preserved here rather than rebranded, and it is
**not** this page's submission. Its measured calibration finding — that any holdout scored against the
public catalogue ranks candidates backwards — is the reason this release also refuses to call its own
proxy a forecast.

### Download the new GeoTIFF — H45

**[Download H45 — the new submission GeoTIFF](docs/downloads/gemsdoe40-h45-eulerdepthreadcluster-20261006-f28e5cff6826-zeros.tif)** · [Live project site](https://buffedlizard55-lab.github.io/GEMSDOE40/docs/index.html) · [H45 executive summary / submission guide](docs/executive-summary-h45.html)

**Unique submission name:** `gemsdoe40-h45-eulerdepthreadcluster-20261006-f28e5cff6826`
**SHA-256 (zeros variant):** `f28e5cff682662c18d6b15f5d5501e13096b0ab00a76e1640dc93dc5e6e6fc8f`

**RESEARCH RELEASE — NOT A PROVEN IMPROVEMENT.** All format checks pass and the file is
measured unique (largest Jaccard 0.0161 against 343 cached prior artifacts, no exact
duplicate). What does *not* exist is a holdout that predicts the organiser's score — see
[the calibration finding](#the-h45-calibration-finding) below. Do not spend a weekly
submission slot on H45 without reading that section.

Upload the `-zeros.tif` variant: every cell is finite and in [0, 1], which is what the
submission form's validator requires. The `-nan.tif` twin is research-only.

### What H45 is

Euler deconvolution depth-clustering over the **magnetic and gravity** layers, exactly as the
brief specifies — *not* gradient thresholding:

1. **Solve.** Reid, Allsop, Granser, Millett & Somerton (1990) Euler deconvolution on
   `tmi` (magnetic, SI 0 and 1) and `iso_grav_anom` (gravity, SI 0) at windows 9/15/25.
   Structural indices come from Reid & Thurston (2014) Table 1. The solver is verified
   against synthetic truth: `tests/test_euler_h45.py`, **10 tests passing**.
2. **Cloud.** 27,533 accepted depth-labelled solutions, each carrying a depth, a conditional
   depth standard error and an Euler-equation misfit.
3. **Weight.** `exp(-z/900)` for shallowness × misfit × depth precision × local depth MAD
   (mutual consistency) × cross-scale corroboration.
4. **Field.** Bilinear mass-conserving splat then a separable Gaussian KDE (σ = 2 px).
5. **Emit.** ~300 m dot lattice (measured median nearest-neighbour 4.0 px), 43,038 dots,
   with the 200 m known-fault ring excluded and the covariate marginals calibrated by
   iterative proportional fitting.

Pipeline entry points: `scripts/run_h45_euler.py` (deconvolution) →
`scripts/build_h45_candidate.py` (emission) → `scripts/audit_h45.py` (uniqueness gate).

### The H45 calibration finding

This is the session's most important output. **Any holdout scored against the public
catalogue ranks candidates backwards.**

| Owner submission | Reported public score | Catalogue lift | Mass inside 200 m ring | Blocked holdout DTI |
|---|---|---|---|---|
| h25-dotted | 0.2600 | 2.179 | 0.1155 | 0.03176 |
| h32-euler  | 0.2649 | 1.810 | 0.0742 | 0.02075 |
| h27-4      | 0.2708 | 1.339 | 0.0299 | 0.01028 |
| h33-b2     | 0.2778 | 0.697 | 0.0000 | 0.00154 |

Both proxy columns are monotone in the *wrong* direction, and this replicates the earlier
32-submission Spearman ρ = −0.79. Predicting the public catalogue is not merely
uninformative — it is harmful. Full numbers:
[`docs/data/h45-calibration-model-20261006.json`](docs/data/h45-calibration-model-20261006.json).

Two effects *are* measurable and both are exploited:

* **The 200 m ring is worthless.** Known faults are masked from scoring, so ring dots earn
  no true positive but still pay the false-positive penalty. Independently confirmed by the
  owner's own record: pruning 2,545 nearest-catalogue dots moved 0.2708 → 0.2778.
* **Dots, not solid lines.** A truth pixel takes `max_x p(x) k(d)`, so neighbours inside one
  300 m kernel do not add. Credit per dot is `s (1 − s/12)` for spacing `s` px, plateauing at
  `s = 6`. The repository's own ledger is an A/B test: `d2.8` → 0.2600 vs `d1.5` → 0.2477.

A forward model fitted to those four scores, `DTI = T / (0.2 M + 0.8 G)` with `G ≈ 5,667`
truth pixels and credit-per-unit-mass `0.0890`, reproduces all four to within 2%. Its
implication is uncomfortable: **all four have the same placement efficiency**, and their
score gaps come almost entirely from total mass and ring mass. H45 sits at or above the
incumbent operating point *provided* its within-shell placement is no worse — and that
cannot be verified without spending a slot.

### Measured Euler results (negative result, kept on the record)

| Family | Accepted solutions | Median depth | Blocked lift vs catalogue |
|---|---|---|---|
| gravity `iso_grav_anom` · SI 0 · w15 | 6,250 | 1,868 m | 0.891 |
| gravity `iso_grav_anom` · SI 0 · w25 | 5,922 | 2,066 m | 1.053 |
| gravity · SI 0 · w15 · one vertical derivative | 0 | — | no window passed QC |
| gravity · SI −1 (corrected finite contact) | 68 | 12 m | degenerate, excluded |
| magnetic `tmi` · SI 0 · w9 | 6,574 | 491 m | 1.002 |
| magnetic `tmi` · SI 0 · w15 | 2,939 | 690 m | 0.832 |
| magnetic `tmi` · SI 1 (thin sheet edge) · w15 | 5,848 | 1,024 m | 0.992 |

Every family sits at 0.83–1.10, i.e. no measurable catalogue-localisation skill. Given the
anti-correlation above this rules the Euler cloud *in* as a candidate generator but provides
no evidence it beats the incumbent within the catalogue shell. Full table:
[`docs/data/h45-analysis-20261006.json`](docs/data/h45-analysis-20261006.json).

Two structural-index corrections were made against the literature and are worth keeping:
Reid & Thurston (2014) correct the gravity index for a finite contact to **−1** (not the
1990 value), while warning it needs a more generalised formulation — measured here, it is
degenerate, giving 68 solutions at ~12 m depth, so it is excluded. And `(z0²−ρ²)/(ρ²+z0²)²`
is a horizontal line of **dipoles**, whose magnetic index is 2 and not 1.

### Earlier arms (preserved, not rebranded)

* **H4 · contact-offset depth-KDE** — [research TIFF](docs/downloads/gemsdoe40-h4-contact-offset-depthkde-20261006-fab9f6619c02.tif) · [guide](docs/executive-summary.html). **HOLD — DO NOT SUBMIT.** 0/16 truth-bearing block wins on the blocked SGMC proxy.
* **H7 · RTP Euler + gravity-gradient context** — [research TIFF](docs/downloads/gemsdoe40-h7-rtp-euler-gravity-context-3d-kde-20261006-998f660f.tif). **HOLD.** Its 0.000265 uses a different proxy version; do not compare it with H4.
* **H40 · depth-cluster / dotted emission** — [NaN twin](docs/downloads/gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-nan.tif) · [zeros twin](docs/downloads/gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-zeros.tif). **HOLD.** Promotion estimates withdrawn by that session.
* **H4-A · contact network + terrain/QFFD** — [closed negative-result report](docs/reports/h4a-negative-result-20261006.md).

### Commands

```bash
. .venv/bin/activate
PYTHONPATH=src python -m pytest tests/test_euler_h45.py -q     # 10 passed
python scripts/run_h45_euler.py  --out work/h45                # deconvolution -> clouds.npz
python scripts/evaluate_h45.py                                 # family skill measurement
python scripts/build_h45_candidate.py --out work/h45 --tag h45 # emission
python scripts/audit_h45.py work/h45/h45-euler-depthcluster-zeros.tif     --out work/h45/h45-audit-zeros.json                        # uniqueness + format gate
```

### Remaining limitations and next-session work

1. **The holdout is the blocker.** Until a proxy that correlates positively with the
   organiser's score exists, candidate ranking is guesswork. The most promising untried
   idea is to validate against *DEM scarps outside the catalogue* rather than against the
   catalogue itself — the four reference submissions all place only 2.4% of their mass in
   the top-5% relief class versus a 5.2% area baseline, so that axis is currently unexploited.
2. **The forward model is fitted from four owner-reported public scores** on an assumption
   (`T + F ≈ M`) that holds only when distinct faults' kernels do not overlap. Refit it as
   soon as another scored submission exists.
3. **Gravity SI = −1 needs the generalised formulation** that Reid & Thurston (2014) say it
   requires; the standard form is degenerate on this data.
4. **Gravity at one vertical derivative produced zero accepted windows** — the QC gates
   reject everything. Worth one diagnostic run with relaxed gates.
5. **`work/bench/h28-edge.tif` is a 162-byte error page**, not a raster. Refetch before any
   comparison that uses it.
6. **Chunks 2–6 of Reid et al. (1990) and 1–3 of Reid & Thurston (2014) remain unread.**


## Reproduce on CPU

Python 3.11 was used. No GPU, trained neural network, password or manual data placement is required for this Euler pipeline. It is **not** the previously claimed ready-to-train U-Net pipeline.

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-research.txt
python -m pip install -e .
bash scripts/download_competition_data.sh          # hash-pinned public mirror of the four inputs
python scripts/prepare_data.py                     # writes evidence/data_pins.json
python scripts/fetch_prior_cache.py --inventory docs/data/prior-inventory-20261006.json --cache data/prior
python scripts/calibrate_live_instrument.py        # anchor measurements (live scores are owner-reported)
python scripts/fit_live_transfer.py                # four transfer instruments, leave-one-out tested
python scripts/run_h41_candidate.py --rank gated --flank 1 --publish --cache work/h41_full.pkl --slug h41-msst-gated27
python scripts/audit_h41_candidate.py              # independent audit of the published bytes
python scripts/write_h41_report.py                 # renders docs/reports/h41-results-20261006.md
python scripts/build_h41_site.py                   # renders docs/*.html from the receipts
python -m pytest -q
python scripts/build_site.py --check
```

The H41 generation step is deterministic: two independent runs of the frozen configuration
produced the identical pixel digest `923c57ab…` and the identical file SHA-256 `1a64839d…`.

Generation first writes to ignored `work/`, reads **no prior prediction or holdout raster**, and does not use a gradient threshold. Only the separate audit runner can publish a research file after exact-format and full raw-output novelty checks. It evaluates a fixed method without a holdout parameter search and never uploads to the competition. For a **new** experiment, refresh the corpus with `scripts/refresh_prior_inventory.py`, then `scripts/audit_prior_history.py --refresh-heads`, review the exclusions, and audit every byte before claiming novelty; the reproducibility command deliberately uses the frozen dated corpus.

Large inputs/prior caches stay out of Git. Only the deliverable TIFF, compressed cloud, compact evidence and actual-data previews are committed. See [three-pass review](docs/reports/review-20261006.md) and [reproduction receipt](docs/data/h4-reproduction.json). The scheduled Pages workflow refreshes **USGS/GDR context only**, exposes stale/error states, and never changes model values or submission eligibility. Automated DrivenData monitoring is not enabled under its [Terms of Use](https://www.drivendata.org/termsofuse/).

## What remains blocked / next session

- **The promotion gate is not cleared.** H41 beats the owner's two best scored files on every instrument this repository can recompute, but H40-E (prior, never scored, sibling repository) leads the blocked proxy mean (0.1395 vs 0.1366) and the best transfer instrument only reaches leave-one-out Spearman 0.705 against a pre-registered 0.80 bar. The next session must either (a) find a genuinely independent validation signal — the two candidates that could move this are a blinded 16-anchor re-fit that survives the bar, or a public expert-mapped fault set that is *not* used in training — or (b) accept the risk explicitly and spend a slot, recording the returned score against the file SHA-256.
- **H40-E is the honest comparator to beat.** It is a supervised detector trained end-to-end on the same surrogate. If a future session can show its proxy advantage is tautological (for example by measuring its calibratable transfer against a held-out fault source), that changes the promotion picture; nothing here establishes it either way.

- **No independent new-fault holdout or hidden organizer labels.** Register a feasible, source-independent whole-fault-system validation protocol before another model search. Do not convert the circular 0.8359 proxy score into a leaderboard forecast.
- **No DrivenData account session / enrollment here.** No authenticated official download or competition submission was made. No credentials are requested. Hash-pinned public mirrors solved data placement, not source authentication.
- **Observation-surface uncertainty.** Audit official GeoDAWN drape/altitude metadata and native gravity resolution before a ground-depth or dip-to-surface projection.
- **Fault versus contact versus resource.** Euler clusters can be unfaulted lithologic boundaries. Finite-step modeling and independent structural/thermal evidence are needed. The competition target is faults, not vents or commercial reservoirs.
- **Score improvement is unproved.** H4 is a valid, genuinely new *research submission-format candidate*, not a winning candidate. The next experiment must beat a credible independent comparator before spending a slot.

## Concurrent mainline work was preserved, not overwritten

`main` moved four times while this release was in review (H8/H13 lineament arms, the browser-QA
hardening, and the H45 depth-read arm). This branch was rebuilt as an overlay on the latest `main`:
every sibling artifact — receipts, downloads, modules, tests and guides — is intact, and the diff
adds this release's files plus the three files that must change to point the site and the README at
the current downloadable file (`docs/index.html`, `docs/sources.html`, `README.md`).

Two consequences are recorded rather than hidden:

* The H45 session's own home page became a sibling page (`docs/executive-summary-h45.html`), still
  linked from every page's navigation and from the sibling-release panel on the overview.
* Only the current release's builders are wired into CI (`scripts/build_h41_site.py --check` and
  `scripts/write_h41_report.py --check`). The retained H4/H8/H13 builders regenerate these same
  five pages and the downloads manifest, so running them now would clobber the current release;
  their receipts remain in the repository and their releases remain downloadable.

The 320-pixel overflow that turned the first CI run of this release red was fixed in
`docs/assets/contact.css` (intrinsic sizing of the numeric checkstrip, which the shared stylesheet
did not constrain), the browser QA now names the offending element in its failure annotation, and
`tests/test_site.py` pins both the rule and a static "no unbreakable token outside a wrapping
context" guard so the class of defect cannot return unnoticed.

## History and operating record

* **Session 2 (this record).** H41 released: Euler SI = 0 depth clusters turned into a unique,
  all-finite, portal-ready GeoTIFF; four live-score transfer instruments fitted and
  leave-one-out tested (best 0.705 — below the bar, published as such); the 0.2778 question
  answered from measurements rather than narrative (mass effect, arrangement effect, pruning
  identity, surrogate anti-correlation); five ranked hypotheses registered; the "[0, 1]" portal
  error traced to two verified causes; site rebuilt from receipts with the download in the first
  viewport; PR raised and merged; remaining work and blockers listed above.

The previous H1/H2 feasibility stops and H2-B negative result remain intact: [append-only validation](docs/research/validation.md), [previous README](docs/archive/README-h2b-20261005.md), [historical H2-B page](docs/index-h2b-20261005.html). Archived mainline outputs are educational references, never this session's candidate. The current site is rebuilt from current JSON receipts; past scores are not silently overwritten.

---

## Full retained user briefs

### Second sitting — 6 October 2026 (H41 release)

<details>
<summary>Read the second brief verbatim (unique TIF, 0.2778 analysis, hypotheses, format bug, site, PR)</summary>

<!-- BEGIN USER BRIEF 20261006B -->
Review the repo. I see that the previous work has been done. But I still don't see the unique tif submission file that I need submitted to the competition. THIS IS THE HIGHEST URGENCY AND MUST BE FOLLOWED!

Also I need the answers to the following questions:

1. Why did 0.2778 score the best on the site and can we get a submission to score even higher than 0.2778? Currently the best score on the leaderboard is 0.3195. I need a submission that is unique and can score higher than 0.2778, ideally close to or above 0.3195.
2. Generate 3-5 candidate geological hypotheses not yet tried, each naming the specific layer(s), the physical signature targeted, why it catches a fault MISSING from the USGS/INGENIOUS catalogue rather than one already in it, and how it differs from anything already in the repo. Rank by expected DTI improvement and implementation cost. Validate the top candidate on our spatially-blocked holdout set before spending a weekly submission slot — do not spend a submission slot on an idea that hasn't beaten the current holdout best. If a candidate needs new external data, name the specific free/official source and confirm it is obtainable before calling the idea viable.
3. Fix the submission-format bug. The downloaded doc failed the form with "Predicted values must be in range [0, 1]" — outputs must be finite and in [0,1] inside the footprint. Provide a unique submission name + short note (e.g. "clustering with k=25").
4. Create a clean, user-friendly GitHub Pages site with the downloadable TIF obvious at the very top, plus an executive-summary subpage explaining exactly how to submit; keep an up-to-date sources/feed section so manual checking is not needed.
5. Store the full prompt/brief in the repo README and re-read it every session; record the Core Values Maximize P(Win) and Own the Outcome as the focal decision framework.
6. Work autonomously, no manual input, multiple passes (Pass 1 implement; Pass 2 review bugs/edge cases; Pass 3 re-check all requirements), line-by-line verification from official trusted sources with links for manual review, flag irregularities, NO HALLUCINATIONS.
7. Create a pull request and merge the pull request onto the main. List the remaining work and any blockers for the next session.
8. Data placement: `bash scripts/download_competition_data.sh` then `python scripts/prepare_data.py`. The repo's Euler pipeline is CPU-runnable; the previously claimed ready-to-train U-Net is not.

Standing constraints: never request passwords/tokens/2FA; do not schedule scraping of DrivenData (its terms of use forbid automated monitoring); fail closed — a candidate that cannot be shown unique, format-valid and better than the current holdout best is not promoted and no weekly slot is spent on it; candidates land in the ignored `work/` directory first and only audited bytes are published; a single submitted file must serve both prize rounds (initial and final).
<!-- END USER BRIEF 20261006B -->

</details>

### First sitting — 6 October 2026 (H4 release, retained)

The complete brief is retained below as the project starting point. Statements in the quoted brief (including historical leaderboard values, GPU requirements and prior readiness claims) are **user-provided context**, not newly verified facts; corrections and limits above take precedence as the measured project record. [Standalone copy](docs/user-prompt-20261006.md).

<details>
<summary>Read the full original task, competition links, score list, constraints and core values</summary>

<!-- BEGIN USER BRIEF 20261006 -->
Review the repo. 

THE FOLLOWING IS THE HIGHEST URGENCY AND MUST BE FOLLOWED!

MUST GENERATE A UNIQUE TIF SUBMISSION FOR THE COMPETITION.  DO NOT COPY A PREVIOUS SUBMISSION UNLESS IT'S FOR LEARNING AND EDUCATION.  BUT WE MUST GENERATE A UNIQUE TIF SUBMISSION.

There should be an easy to download submission tif file as described by the prompt.  Read the entire prompt.

Euler deconvolution depth-clustering. Build this from potential-field depth estimation, not gradient thresholding: run Euler deconvolution (Reid, Allsop, Granser, Millett, and Somerton, Geophysics, 1990) across the magnetic and gravity layers with the structural index for a fault-like contact, producing a cloud of depth-labeled solution points rather than a single edge map. Convert that cloud into a continuous raster via kernel-density estimation of solution density per pixel, weighted so tight clusters of shallow, mutually-consistent solutions score higher than scattered or deep ones, since a real near-surface fault produces the former and noise produces the latter. Normalize to [0,1] and write a single-band float32 GeoTIFF in EPSG:32611, 100 m resolution, matching the sample submission's exact shape and geotransform, NaN only outside the valid footprint — then, before presenting it for download, hash and correlate it against every prior submission's raw output and refuse to call it new if it's a near-duplicate, since depth-clustering should produce a visibly different spatial pattern than any gradient or curvature candidate already made.

The following sites should serve as a starting point for understanding how to generate TIF submissions.  These websites are researched, and tested and have generated TIF submissions.  But we need to generate high scoring submissions.

Here are the results from submissions into the competition, separated by ....:

[https://buffedlizard55-lab.github.io/GEMSDOE/docs/index.html](https://buffedlizard55-lab.github.io/GEMSDOE/docs/index.html)

gems-submission-20260925T001403Z-7f00890a: 0.1563

....

[https://buffedlizard55-lab.github.io/6GEMSDOE/](https://buffedlizard55-lab.github.io/6GEMSDOE/)

gems6_hgb88-topk03_33cec71ff0: 0.0286

....

[https://buffedlizard55-lab.github.io/GEMSDOE3/docs/index.html](https://buffedlizard55-lab.github.io/GEMSDOE3/docs/index.html)

pindrop-v4-nodes-20260925T152420Z-f347b70daa: 0.1193

pindrop-v4-discovery-20260925T152423Z-37f9d5b855: 0.0830

pindrop-v4-ridge-20260925T152422Z-4e03fc9705: 0.1152

....

[https://buffedlizard55-lab.github.io/GEMSDOE2/docs/index.html](https://buffedlizard55-lab.github.io/GEMSDOE2/docs/index.html)

gemsdoe2-dual-family-union-20260925T160406Z-f68e590f: 0.1560

....

[https://buffedlizard55-lab.github.io/GEMSDOE4/](https://buffedlizard55-lab.github.io/GEMSDOE4/)

gems-submission-20260926T163915Z-237f0063: 0.0343

....

[https://buffedlizard55-lab.github.io/5GEMSDOE/docs/index.html](https://buffedlizard55-lab.github.io/5GEMSDOE/docs/index.html)

gems-submission-20260926T175114Z-7f00890a: 0.1563

....

[https://buffedlizard55-lab.github.io/7GEMSDOE/](https://buffedlizard55-lab.github.io/7GEMSDOE/)

lidarscarp-ridge-top2pct-36c3a3f341c8: 0.1461

....

[https://buffedlizard55-lab.github.io/8GEMSDOE/](https://buffedlizard55-lab.github.io/8GEMSDOE/)

Hedge-v2_submission: 0.1563

....

[https://buffedlizard55-lab.github.io/GEMSDOE9/docs/index.html](https://buffedlizard55-lab.github.io/GEMSDOE9/docs/index.html)

2314b599: 0.0107

....

[https://buffedlizard55-lab.github.io/11GEMSDOE/docs/index.html](https://buffedlizard55-lab.github.io/11GEMSDOE/docs/index.html)

gems-structural-area06-v1: 0.0202

....

[https://buffedlizard55-lab.github.io/12GEMSDOE/docs/index.html](https://buffedlizard55-lab.github.io/12GEMSDOE/docs/index.html)

r7-nms3-dem10-scarp_0c9199f14e62:0.1294

r7-nms3-dem10-scarp_0c9199f14e62_allfinite:0.1294

....

[https://buffedlizard55-lab.github.io/15GEMSDOE/docs/index.html](https://buffedlizard55-lab.github.io/15GEMSDOE/docs/index.html)

gems-tso1-20260929T005627Z-conj_alteration_mag: 0.0782

....

[https://buffedlizard55-lab.github.io/14GEMSDOE/docs/index.html](https://buffedlizard55-lab.github.io/14GEMSDOE/docs/index.html)

GEMS_r5-geom-horse-ensemble_20260929T154852Z_ccbe1de0_site_e96e942f: 0.0020

....

[https://buffedlizard55-lab.github.io/17GEMSDOE/](https://buffedlizard55-lab.github.io/17GEMSDOE/)

17GEMSDOE_F-ensemble-2pct_20260930T050626Z:0.0187

....

[https://buffedlizard55-lab.github.io/18GEMSDOE/](https://buffedlizard55-lab.github.io/18GEMSDOE/)

H19-C_20260930T212401Z_c11e495e: 0.0297

....

[https://buffedlizard55-lab.github.io/19GEMSDOE/docs/index.html](https://buffedlizard55-lab.github.io/19GEMSDOE/docs/index.html)

h19-4-multiline-corroborated-openness-thermal-pop-20260930-691e4dfa-nan: 0.1894

h19-5-powerlaw-budget-multiline-corroborated-20260930-e27054cf-nan: 0.1922

....

[https://buffedlizard55-lab.github.io/GEMSDOE10/](https://buffedlizard55-lab.github.io/GEMSDOE10/)

h16-continuation-20260927T065521077735Z-3431b83c7c: 0.0461

h20-dem10-scarp-thin-20260927T155223039488Z-ffc91a1686: 0.0921

H25-ctx-ridge-20260927T232947704150Z-6452ae1d00: 0.1280

h28-dotted-ridge-20260928T020256236880Z-6452ae1d00: 0.1839

....

[https://buffedlizard55-lab.github.io/13GEMSDOE/](https://buffedlizard55-lab.github.io/13GEMSDOE/)

20261001_r13-lattice-s5_v2_nan-outside:0.0904

....

[https://buffedlizard55-lab.github.io/16GEMSDOE/docs/index.html](https://buffedlizard55-lab.github.io/16GEMSDOE/docs/index.html)

h16-1-topo-geophys-baseline-ridges-20260930-df20f65e-nan: 0.1855

h18-3a-topo-geophys-x-complexity-prior-20260930-c502dfab-nan: 0.0976

h18-4-usgs-geologic-map-faults-gap-20260930-aef8f42c-nan: 0.0360

....

[https://buffedlizard55-lab.github.io/GEMSDOE21/](https://buffedlizard55-lab.github.io/GEMSDOE21/)

h19-4-reference-20260930-691e4dfa: 0.1894

....

[https://buffedlizard55-lab.github.io/20GEMSDOE/docs/index.html](https://buffedlizard55-lab.github.io/20GEMSDOE/docs/index.html)

h20-1-sarnnpu-powerlaw-pi0363-tilt-wingcrack-20260930-be0e8f6b-nan: 0.1890

h20-5-continuous-pu-proxy-unverified-20260930-824ce73a-nan: 0.1859

....

[https://buffedlizard55-lab.github.io/GEMSDOE22/docs/index.html](https://buffedlizard55-lab.github.io/GEMSDOE22/docs/index.html)

h23-a-dti-optimal-emission-6pct-20261002-e2ec4b49-nan: 0.1002

h23-b-dti-optimal-emission-10pct-20261002-86176698-nan: 0.0748

....

[https://buffedlizard55-lab.github.io/GEMSDOE23/](https://buffedlizard55-lab.github.io/GEMSDOE23/)

h30-arrangement-matched-habitat-20261002-0d4e02e8-nan: 0.1352

....

[https://buffedlizard55-lab.github.io/GEMSDOE24/](https://buffedlizard55-lab.github.io/GEMSDOE24/)

h25-1-dotted-h19-5-d1-5-20261002-989f59505db1-nan: 0.2477

....

[https://buffedlizard55-lab.github.io/GEMSDOE25/](https://buffedlizard55-lab.github.io/GEMSDOE25/)

dotted-h19-5-d2-8-20261002-e56ea318af89-nan: 0.2600

....

[https://buffedlizard55-lab.github.io/GEMSDOE26/](https://buffedlizard55-lab.github.io/GEMSDOE26/)

dilcond-oof-v1-20261003-47629f496133-nan: 0.1223

....

[https://buffedlizard55-lab.github.io/GEMSDOE27/](https://buffedlizard55-lab.github.io/GEMSDOE27/)

topo-gap-closure-t-v2-on-d1-5-20261002-5512495c6bd1-nan: 0.2449

....

[https://buffedlizard55-lab.github.io/GEMSDOE28/](https://buffedlizard55-lab.github.io/GEMSDOE28/)

h27-4-r1-solo-d2-8-20261003-8acb75e1f2cc-nan: 0.2708

h32-1-prethin-tip-euler-d2-8-20261003-31e35eee884e-nan: 0.2649

h36-1-rung30-blind-r1-20261003-b531dae0a36f-nan:

h38-1-hf-euler-r30-r1-20261003-56a9f473edc7-nan:

....

[https://buffedlizard55-lab.github.io/GEMSDOE29/docs/index.html](https://buffedlizard55-lab.github.io/GEMSDOE29/docs/index.html)

efd28-repro-20261003-1cc7dc534d51-nan: 0.2600

repo-c0-habitat-emission-20261003-a4d439b07426-nan: 0.0041

sgmc-off-catalogue-44k-20261003-c8dcd780e3fd-nan: 

wormrank-d28-20261003-59dcaf6dd11d-zeros:

wormsurv-filter-20261003-921f10960d6e-zeros:

xfit-c0-habitat-20261003-ca879db0089a-zeros:

xfit-h41-union-qfaults-20261003-9edb34b99e3a-zeros:

....

[https://buffedlizard55-lab.github.io/GEMSDOE30/](https://buffedlizard55-lab.github.io/GEMSDOE30/)

d28-poisson300m-offcat-44090-20261003T233156Z-91eae1ca: 0.2600

....

[https://buffedlizard55-lab.github.io/GEMSDOE31/docs/](https://buffedlizard55-lab.github.io/GEMSDOE31/docs/)

h27-4-solo-d28-20261004-8acb75e1-nan:0.2708

....

[https://buffedlizard55-lab.github.io/GEMSDOE32/docs/index.html](https://buffedlizard55-lab.github.io/GEMSDOE32/docs/index.html)

h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros: 0.2778

....

[https://buffedlizard55-lab.github.io/GEMSDOE33/](https://buffedlizard55-lab.github.io/GEMSDOE33/)

h33d-analog-tip-stepover-r30-20261004-cb490425926e: 0.2632

....

[https://buffedlizard55-lab.github.io/GEMSDOE34/docs/index.html](https://buffedlizard55-lab.github.io/GEMSDOE34/docs/index.html)

h34-scatter-q50-arr-matched-20261004T223317Z: 0.0778

....

[https://buffedlizard55-lab.github.io/GEMSDOE35/docs/index.html](https://buffedlizard55-lab.github.io/GEMSDOE35/docs/index.html)

h35-06-aaa86efb25-20261004T225420098147Z-candidate: 0.0418

....

[https://buffedlizard55-lab.github.io/GEMSDOE36/docs/](https://buffedlizard55-lab.github.io/GEMSDOE36/docs/)

anderson-geothermal-pinn-38854-20261004T230000Z-9b9ea4e6-zeros: 

....

[https://buffedlizard55-lab.github.io/GEMSDOE37/](https://buffedlizard55-lab.github.io/GEMSDOE37/)

h6-physics-dotted-80k-20261005T055000Z-0bef9211631c:

....

[https://buffedlizard55-lab.github.io/GEMSDOE38/docs/index.html](https://buffedlizard55-lab.github.io/GEMSDOE38/docs/index.html)

D-step-3p0-07pct-tipProt-20261005-ecfbf59e2b48-zero:

....

[https://buffedlizard55-lab.github.io/GEMSDOE39/](https://buffedlizard55-lab.github.io/GEMSDOE39/)

h40-e-disc-h40e-30k-zeros:

....

40GEMSDOE

:

....

41GEMSDOE

:

....

42GEMSDOE

:

....

43GEMSDOE

:

....

44GEMSDOE

:

....

WE NEED TO STUDY, ANALYZE, AND UNDERSTAND THE HIGHEST SCORE FROM THE GEMDOE SITE WHERE THE SUBMISSION TIF IS DOWNLOADED FROM WHICH IS THE FOLLOWING:

[https://buffedlizard55-lab.github.io/GEMSDOE32/docs/index.html](https://buffedlizard55-lab.github.io/GEMSDOE32/docs/index.html)

h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros: 0.2778

Why and how did this get the highest score and are we able to generate a submission that scores higher than 0.2778?

Answer the question using Phd level experience, knowledge, and judgement. Then use the answer to generate a unique TIF submission into the competition.  Must be unique submission unlike any within the GEMSDOE sites above.  Verify working line by line no hallucinations.

The following is the leaderboard for the competition:

[https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/](https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/)

We need to quickly look at the results and results from the GEMSDOE websites above.

Before implementing, generate 3–5 candidate geological hypotheses we haven't tried yet, each naming: the specific layer(s) involved, the physical signature being targeted (e.g., an edge-detection or curvature transform), why it should catch a fault missing from the USGS/INGENIOUS catalogue rather than one already in it, and how it differs from anything already implemented in this repo. Rank them by expected DTI improvement and implementation cost. Validate the top candidate on our spatially-blocked holdout set before touching a weekly submission slot — do not spend a submission slot on an idea that hasn't beaten the current holdout best. If a candidate can't be validated without new external data, name the specific free, official source needed and check it's obtainable before proposing the idea as viable.

Work line by line verifying from official verified trusted sources, provide links for manual review.  There should be no manual input, work on your own to complete tasks.  Flag any irregularities for review.  No hallucinations.                      

Verify no hallucinations.    

The goal of this project is to get a full list that follow our requirements.  No hallucinations.  Verify line by line.

We have a good understanding of how our hypothesis, methodology, calculations, analysis are done so we should be able to figure out a way to score higher on the leaderboard using previous results and scoring that we have across the sites listed above.  We need to come up with distinct and unique strategies to score higher in this competition leaderboard.  We need to start doing heavy and deep research into the part of the project that matters the most, which is the scientific discovery of geothermal vents.  We should store all of our information and knowledge that we can gather from official verified sources.  This will serve as a starting point for other projects as well.  We need to think outside the box but still be grounded in proper scientific research, we are ultimately aiming for a top prize that many others are competing for.  So it's important to be contrarian but be smart about it.  We need to find sources of data that others are over looking or areas of the project when it comes to geothermal vents.  We need to do deep research and critical thinking and come up with new hypothesis to test.

0.3195	is the highest score right now so we need to design a new strategy, research, testing, analyzing, and generating submission system than the current website.  It should be unique, take unique approaches to generating a submission that can score higher than 0.3195.  

Put this prompt into the repo readme and read it everytime we work on the project as a starting point to make sure we are building what we are aiming for and have a strong base to continue building and improving on making something useful for everyday use.  It should solve the problem of having to manually check everything ourselves and having an up to date current feed.

Review the repo. 

The following is taken from the Arena AI team and I think it makes a good point on building a successful project, so let's keep the Core Values and Own the Outcome as a focal point when building, developing, researching, suggesting upgrades, and implementing the work.

Our Core Values

Maximize P(Win)

“Maximize the Probability of Winning”: our decision making framework. In every decision, we weigh tradeoffs, assess risk, and choose the path that maximizes the probability that Arena succeeds. We set aside our emotions and make tough decisions in order to maximize P(Win). “Maximize P(Win)” frees us from constraints and clarifies that we must put Arena first.

Own the Outcome

We own results end to end — not just our individual slice of the work. When problems arise and we have the means to act, we do so without waiting for permission or assignment. We treat failure and success as signals and use them to improve. At Arena, we stay accountable to the final outcome.

Work line by line verifying from official verified trusted sources, provide links for manual review.  There should be no manual input, work on your own to complete tasks.  Flag any irregularities for review.  No hallucinations.                      

Verify no hallucinations.    

The goal of this project is to get a full list that follow our requirements.  No hallucinations.  Verify line by line.

We need to focus on being able to generate a submission into the competition.  

The site should be able to generate a TIF file that is required for submission.  It should be as easy as download to click a File to submit into the competition.  This needs to be in the executive summary or the very beginning of the site.  it should be obvious when you visit the site.

I tried to submit the document that i downloaded from the site but it returned this error on the submission form:

"Predicted values must be in range [0, 1]"

Also we need to give it a unique name and A short comment to help you or your team tell submissions apart later e.g. clustering with k=25

Here is the submission page when i click submit file

New submission

File to submitNo file chosen

You can submit a single-band GeoTIFF (.tif) file, or a .zip file containing a single GeoTIFF, with your predictions. It must match the submission format's CRS, shape, and geotransform. You may wish to review the competition rules first.

Note (optional)

A short comment to help you or your team tell submissions apart later e.g. clustering with k=25

Create a executive summary subpage that explains exactly how to make a submission into the contest.

Work on the next steps from the previous sessions first.

The goal of this project is to place top of the leaderboard in this competition.  The following is the competition:

[https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/)

We need to create a project that can compete and place top of the leaderboard.  We need to understand the problem, collect all the data and organize it into a clean easily auditable table with official verified links for manual verification.  

This is the guidelines we need to follow.[https://www.drivendata.org/competitions/306/competition-doe-gems/](https://www.drivendata.org/competitions/306/competition-doe-gems/)

Get familiar with the problem through the overview and problem description,[https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/). You might also want to reference additional resources available on the about page,[https://www.drivendata.org/competitions/306/competition-doe-gems/page/968/](https://www.drivendata.org/competitions/306/competition-doe-gems/page/968/).

Download the data from the data,[https://www.drivendata.org/competitions/306/competition-doe-gems/data/](https://www.drivendata.org/competitions/306/competition-doe-gems/data/), tab.  

Create and train your own model. This reference solution,[https://github.com/drivendataorg/gems-prize-reference-solution](https://github.com/drivendataorg/gems-prize-reference-solution) implements a simple approach.

Use your model to generate predictions that match the submission format.

Tell me what are you limitations and what you need access to during this project.  We will need to find free publicly available sources and data from official and verified sources if we are to use 3rd party or external data.  

this pdf outlines how submissions must be entered into the competition.  

[https://docs.nlr.gov/docs/fy26osti/96647.pdf](https://docs.nlr.gov/docs/fy26osti/96647.pdf)

You must be able to do your own research, deep research, scientific literature research and organize the knowledge so that we can critically think through the problem and generate a solution through scientific and free publicly available information.  this must be done autonomously and must be constantly reviewed and improved upon.  Provide suggestions and improvements and implement them.

❌ No DrivenData auth → cannot auto-download training_features.tif, labels.tif, sample_submission.tif, 1m_DEM_links.csv from [https://www.drivendata.org/competitions/306/competition-doe-gems/data/](https://www.drivendata.org/competitions/306/competition-doe-gems/data/) (verified redirect to login)

See below for links from the above site.  See attached files for links from the above site.

[https://gdr.openei.org/submissions/1391](https://gdr.openei.org/submissions/1391)

Download competition data from [https://www.drivendata.org/competitions/306/competition-doe-gems/data/](https://www.drivendata.org/competitions/306/competition-doe-gems/data/) (requires login) to data/

See links below for competition data:

[https://www.dropbox.com/scl/fi/aemhtutjgcp6tr3tint94/GEMS_96647.pdf?rlkey=rek210cj2smnmzb8n0sla1vmd&amp;st=wz4kofki&amp;dl=0](https://www.dropbox.com/scl/fi/aemhtutjgcp6tr3tint94/GEMS_96647.pdf?rlkey=rek210cj2smnmzb8n0sla1vmd&st=wz4kofki&dl=0)

[https://www.dropbox.com/scl/fi/6rgvnuady818ol8yqgis4/example_submission.tif?rlkey=kbykilvau066xuogoosbf4cq8&amp;st=8junzdyw&amp;dl=0](https://www.dropbox.com/scl/fi/6rgvnuady818ol8yqgis4/example_submission.tif?rlkey=kbykilvau066xuogoosbf4cq8&st=8junzdyw&dl=0)

[https://www.dropbox.com/scl/fi/t7fyt03qdh9egyme0itwo/existing_faults.tif?rlkey=yiao96uluqdkipf0h5vju71jf&amp;st=rnino7ya&amp;dl=0](https://www.dropbox.com/scl/fi/t7fyt03qdh9egyme0itwo/existing_faults.tif?rlkey=yiao96uluqdkipf0h5vju71jf&st=rnino7ya&dl=0)

[https://www.dropbox.com/scl/fi/3vz9o0wwavi26xaeoxlwr/gems-geodawn-numerical-features.tif?rlkey=je8d8fepqfbst9lnwsq9rkplu&amp;st=zj1lag1r&amp;dl=0](https://www.dropbox.com/scl/fi/3vz9o0wwavi26xaeoxlwr/gems-geodawn-numerical-features.tif?rlkey=je8d8fepqfbst9lnwsq9rkplu&st=zj1lag1r&dl=0)

[https://www.dropbox.com/scl/fi/ig0mban712ns1atphgphe/Digital-elevation-model-links-JSON.pdf?rlkey=zm77f1vbtt2if8hlruymptnu3&amp;st=srhhir10&amp;dl=0](https://www.dropbox.com/scl/fi/ig0mban712ns1atphgphe/Digital-elevation-model-links-JSON.pdf?rlkey=zm77f1vbtt2if8hlruymptnu3&st=srhhir10&dl=0)

Work line by line verifying from official verified trusted sources, provide links for manual review.  There should be no manual input, work on your own to complete tasks.  Flag any irregularities for review.  No hallucinations.                      

Verify no hallucinations.    

The goal of this project is to get a full list that follow our requirements.  No hallucinations.  Verify line by line.

Site creation

Create a github page for this repo that has clean ui, user friendly, simple and easy to use.  It should be organized and clean.  

It should include all relevant information in an easy to read format with official verified links as sources for review.  Work line by line verify everything no hallucinations.

**The single remaining blocker to training is data placement**: run `bash scripts/download_competition_data.sh` on any unrestricted machine into `data/`, then `python scripts/prepare_data.py` — after that the full train→inference→validate pipeline is ready to run (GPU needed for training; metric/losses/validation all verified working here on CPU).

you need to complete the above task by yourself.  Work line by line verifying from official verified trusted sources, provide links for manual review.  There should be no manual input, work on your own to complete tasks.  Flag any irregularities for review.  No hallucinations.                      

Verify no hallucinations.    

The goal of this project is to get a full list that follow our requirements.  No hallucinations.  Verify line by line.

Run this task through multiple passes.

Pass 1: Implement the task completely and verify the result.

Pass 2: Review your work for bugs, missing requirements, incorrect assumptions, and edge cases. Fix everything you find.

Pass 3: Re-check the entire implementation against the original request. Improve accuracy, reliability, completeness, and code quality. Fix any remaining issues.

Do not stop after the first pass. Each pass must build on the previous one. Before finishing, verify that the final result fully satisfies the original request.  Work line by line verify everything no hallucinations.

Go ahead and create a pull request and then merge the pull request onto the main. Make suggestions for what work still needs to be done and any limitations that is in the way of a successful project.  It should be worked on in this next session or the next session.  Work line by line verify everything no hallucinations.

<!-- END USER BRIEF 20261006 -->

</details>
