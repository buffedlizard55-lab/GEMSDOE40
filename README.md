# GEMSDOE40 — SI = 0 contact Euler depth-clustering submission candidate

> **Read this README and the complete retained brief below at the beginning of every session.** See [AGENTS.md](AGENTS.md). Principles: **Maximize P(Win)** and **Own the Outcome**.

## Download the unique submission GeoTIFF

**[Download H8 — the SI = 0 contact Euler depth-clustering candidate](docs/downloads/gemsdoe40-h8-euler-lineament-depthcluster-20261006-785c4f5d5ce1.tif)** · [Metric-optimal twin (identical dots, every value 1.0)](docs/downloads/gemsdoe40-h8-euler-lineament-depthcluster-20261006-785c4f5d5ce1-hard.tif) · [Live project site](https://buffedlizard55-lab.github.io/GEMSDOE40/docs/index.html) · [Executive summary / submission guide](docs/executive-summary.html)

**Status: BUILT, AUDITED, NOT PROMOTED.** The format gate passed and the raw-output novelty gate passed against every cached prior raster, but the repository's promotion rule is *not* satisfied: no local instrument ranks the family's recorded scores at ρ ≥ 0.8 (best measured: mass -0.676, SGMC off-catalogue credit +0.676, catalogue-calibrated LM +0.100). **No leaderboard score is claimed for this file and no weekly submission slot was spent.** Uploading is a human decision, not a recommendation of this repository.

| Item | Measured result |
|---|---|
| File | `gemsdoe40-h8-euler-lineament-depthcluster-20261006-785c4f5d5ce1.tif` (493,872 bytes) |
| Grid | EPSG:32611, 100 m, 3,730 rows × 3,292 columns, transform `(100, 0, 243350, 0, -100, 4508550)` — the exact sample transform |
| Values / footprint | 40,000 positive cells; continuous 0.383 … 1.000 with 50 % of dots at exactly 1.0; all 5,167,373 inside cells finite in [0,1]; all 7,111,787 outside cells NaN |
| TIFF encoding | Single band, float32, DEFLATE, **Predictor = 1 (none)** — the official `sample_submission.tif` setting, never the integer predictor 2 that caused the one observed rejection |
| File SHA-256 | `49d26d8195e343b531cf5e8f9ebd44e2e0a652daca6ad4f71016faabcab33509` |
| Canonical pixel SHA-256 | `785c4f5d5ce1b2a341f27ac26702b77626eaed270b6af48e0127d68378983500` |
| Method | SI = 0 contact Euler (Reid et al. 1990, eq. 2 with the arbitrary offset A) on `rtp` (150 m continuation, windows 9/15/21) and once-differentiated `iso_grav_anom` (500 m continuation, windows 15/21/31); stride 4; 236,401 QC-passing solutions |
| Weighting | Shallow-depth decay (1200 m) × cross-window depth consensus (≥ 3 neighbours in 300 m) × lineament coherence (1 − λ₂/λ₁) × cross-family corroboration (400 m / 600 m, uncorroborated ×0.5) |
| Emission | Anisotropic KDE (σ 1.6 along / 0.55 across, 12 direction bins) → value-ranked non-maximum suppression at 2.8 px under a 40,000-dot budget; catalogue ± 1 px excluded |
| Cloud | 124,382 retained solutions (40,185 magnetic + 84,197 gravity) with row, col, depth, depth SE, residual, condition, window, cluster weight, lineament coherence and corroboration flag |
| Novelty | Compared against the 279 cached prior rasters (278 comparable, 1 unreadable); max |Pearson| 0.0253, max top-mass Jaccard 0.0154, 0 exact duplicates |
| Named surrogates | Catalogue w = 0.0174 (cover 0.0178); SGMC off-catalogue w = 0.0376 (cover 0.0453) — proxies, never scores |
| Solver control | Analytic contact recovered to 0.3–1.0 m laterally and exactly in depth at zero noise (585/912/1,417 accepted at windows 9/15/25); 4 solutions survive at σ = 1; none at σ = 5 |
| Decision | **NOT PROMOTED** — `slot_eligible = false`, `organizer_score = null`, `weekly_submission_used = false` |

**Unique tracking name (for your own record):** `GEMSDOE40-H8-LINEAMENT-785c4f5d5ce1`

**Short note for the submission form (188 characters):**

> Euler deconvolution SI=0 contact depth-clustering over magnetic (rtp) and isostatic gravity (differentiated once); lineament-weighted kernel-density emission; 40,000 dots; id 785c4f5d5ce1.

[Depth-labelled solution cloud](docs/downloads/h8-lineament-solutions.csv.gz) · [Audit receipt](docs/data/h8-lineament-audit.json) · [Generation receipt](docs/data/h8-lineament-generation.json) · [Why 0.2778 leads](docs/research/h8-analysis-20261006.md) · [Pre-registration](docs/research/h8-preregistration-20261006.md)

## What the measurements say

1. **The best recorded score is a deletion, not a discovery.** The 0.2778 file is the 0.2708 base minus 2,545 pixels (6.33 % of its mass), every one of them 141–200 m from the provided catalogue, with **0 pixels added**. Its minimum distance to a mapped fault is 224 m. Pruning catalogue-adjacent mass is the entire gain.
2. **Every local promotion instrument fails, and that is measured, not assumed.** Against the sixteen recorded scores: mass ρ = -0.676, SGMC off-catalogue credit density ρ = +0.676, catalogue-calibrated LM ρ = +0.100; a per-block credit field fitted to fifteen scores and tested out of sample scores negatively. The promotion rule therefore cannot be satisfied with local data, and no slot is spent.
3. **What beating 0.3195 would take.** Because α + β = 1, one more unit of predicted mass moves the metric denominator by exactly 0.2 wherever it lands, so a dot only pays while its kernel credit exceeds 0.2 × DTI (≈ 283 m from a hidden fault at DTI = 0.2778). Inverting the metric on the recorded scores puts the best file at ≈ 0.090 mean credit per dot; reaching 0.3195 at the family's own mass needs ≈ 0.103 — roughly 15 % better placement, not exotic new data.
4. **The solver is validated; the data is the limit.** On an exact analytic contact the frozen SI = 0 solver returns the contact to 0.3–1.0 m laterally and exactly in depth; on the real layers the solution cloud sits at chance against the mapped faults (catalogue w = 0.0174). That is a statement about the potential-field contacts available at 100 m, not about the code.
5. **Ranked next hypotheses** (full register in [docs/hypotheses.html](docs/hypotheses.html)): **H9** tip / step-over / along-strike extensions of mapped systems (highest expected gain per cost; testable today with leave-the-tips-out); **H10** finite-step gravity inversion at SI = −1; **H11** 1 m lidar scarp re-mapping — *flagged irregularity*: the earlier scarp stack is absent from this workspace and the USGS S3 endpoint is TLS-blocked from this sandbox, so it cannot be rebuilt or re-measured here; **H12** microseismicity alignment (cheapest, weakest measured skill).
6. **Retained history stays downloadable and held.** H4 (`gemsdoe40-h4-contact-offset-depthkde-20261006-fab9f6619c02.tif`), H7 (`gemsdoe40-h7-rtp-euler-gravity-context-3d-kde-20261006-998f660f.tif`) and H40 run-2 are published with explicit **HOLD — DO NOT SUBMIT** labels; their byte hashes are pinned in `scripts/build_site.py` and verified by CI, so a silent change cannot pass.

## Reproduce on CPU

Python 3.11, no GPU, no learned weights, no randomness anywhere in generation, no competition login.

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-research.txt
python -m pip install -e .
python scripts/acquire_data.py          # sha256-pinned data bridge, no login wall in this sandbox
OPENBLAS_NUM_THREADS=2 python scripts/run_h8_lineament.py   # writes work/h8lineament/*
OPENBLAS_NUM_THREADS=1 python scripts/audit_h8_lineament.py --publish   # audit + publication gate
python scripts/build_h8_site.py         # rebuild the five site pages from the receipts
python -m pytest -q && python scripts/build_site.py --check
```

## Repository map

`src/gemsdoe40/contact_euler.py` frozen contact solver · `src/gemsdoe40/h8_euler.py` H8 weighting + anisotropic KDE ·
`src/gemsdoe40/h8_solver.py` local window loop · `scripts/run_h8_lineament.py` generation ·
`scripts/audit_h8_lineament.py` format/novelty/surrogate audit and publication gate ·
`scripts/build_h8_site.py` site generator (wrapped by `scripts/build_contact_site.py` for CI compatibility) ·
`docs/research/h8-analysis-20261006.md` metric algebra · `docs/research/h8-preregistration-20261006.md` frozen settings ·
`docs/data/*.json` machine-readable receipts · `data/prior/` byte-verified prior rasters (gitignored).

## Limitations that remain in the way

- **No organizer score exists for any published file here.** Leaderboard rows are owner-reported and are not per-file receipts.
- **No local instrument can validate a submission.** The promotion gate is unmeetable without a truth model that ranks the recorded scores (ρ ≥ 0.8), which is the first task of the next session.
- **The competition data is login-walled.** This sandbox builds from a sha256-pinned bridge; it cannot download or verify DrivenData files itself.
- **1 m lidar is unavailable here.** The 3DEP endpoint is TLS-blocked from this sandbox, and the previously built scarp stack is absent, so H11 is registered but not measurable in this environment.
- **A potential-field contact is not a fault.** It may be a lithologic boundary, and a fault is not automatically a geothermal reservoir.


---

## Full retained user brief — 6 October 2026

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
