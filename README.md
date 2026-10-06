# GEMSDOE40 — contact Euler depth-clustering

> **Read this README and the complete retained brief below at the beginning of every session.** See [AGENTS.md](AGENTS.md). Principles: **Maximize P(Win)** and **Own the Outcome**.

## Download the new GeoTIFF

**[Download H4 — the new submission-format research TIFF](docs/downloads/gemsdoe40-h4-contact-offset-depthkde-20261006-fab9f6619c02.tif)** · [Live project site](https://buffedlizard55-lab.github.io/GEMSDOE40/docs/index.html) · [Executive summary / submission guide](docs/executive-summary.html)

**HOLD — DO NOT SUBMIT.** The file is genuinely distinct and format-valid, but it did **not** beat the blocked holdout best. No weekly slot was used; no organizer score is known. Publishing a valid research download is not a recommendation to upload it.

| Item | Measured result |
|---|---|
| File | `gemsdoe40-h4-contact-offset-depthkde-20261006-fab9f6619c02.tif` |
| Name | `GEMSDOE40-H4-CONTACT-OFFSET-fab9f6619c02` |
| Format | Single-band float32, EPSG:32611, 100 m; 3,730 rows × 3,292 columns |
| Affine transform | `(100, 0, 243350, 0, -100, 4508550)` — exact sample transform |
| Range / footprint | All 5,167,373 inside pixels finite in [0,1]; all 7,111,787 outside pixels NaN |
| File SHA-256 | `ee73ffd76fabbaa1a2e77e17f57947a7db858916d713801e0c3502e49f6acabb` |
| Canonical pixel SHA-256 | `fab9f6619c02595ceec668c87420789035d8103d3c93b91a2e4aa7c2bc905e86` |
| New construction | SI=0 **with contact offset A**, rank-adaptive TMI + first-vertical-gravity-derivative Euler; depth-consistent cross-window KDE |
| Cloud | 46,656 QC-passing solutions; **12,167 magnetic + 8,450 gravity** clustered solutions |
| Continuous output | 865,145 positive cells; 857,785 distinct finite values; confidence mass 29,350.41 |
| Prior-output audit | **343 byte-verified artifacts** across 55 public repositories / 338 pinned branch heads / 1,859 reachable commits; 339 exact-grid, 340 full-size array comparisons, 2 supplementary alignment checks; 190 ZIP paths/versions inspected. One spatial format fixture, a one-pixel demo and a 110-byte header-only non-raster are explicitly distinguished |
| Novelty | No exact/near duplicates; maximum absolute Pearson **0.128319**, top-37,654 Jaccard **0.049527**, containment **0.137676** |
| Local SGMC proxy DTI | H4 **0.006902**; H33-B2 **0.091550**; frozen SGMC-derived incumbent **0.835907** |
| Decision | 0/16 truth-bearing block wins vs both H33 and incumbent; **HOLD — DO NOT SUBMIT** |

**Short note (137 characters, research record only):**

> H4 SI0+A rank-aware TMI+dGdz Euler; shallow cross-window depth-consensus KDE. New raw field; HOLD, proxy gate failed; no organizer score.

[Download all depth-labeled solutions](docs/downloads/h4-euler-solutions.csv.gz) · [Exact format receipt](docs/data/h4-format.json) · [Every prior comparison](docs/data/h4-uniqueness.json) · [Every holdout block](docs/data/h4-validation.json) · [Human-readable result](docs/reports/h4-results-20261006.md)

## What the review found

1. **The requested method was not already correct.** The older H1 contact solve forced A=0, while [Reid et al. (1990), equation 2 / step 3b](https://www.reid-geophys.co.uk/wp-content/uploads/2017/11/Reid-et-al-1990.pdf) explicitly fit an arbitrary contact offset. A full-rank position requirement also discards ideal straight contacts. H4 tests these cases analytically and handles their unresolved along-strike coordinate explicitly.
2. **Gravity needs a caveat.** SI=0 on raw gravity is not the magnetic-contact model. H4 uses SI=0 on dG/dz as a **local top-edge approximation**, not an exact finite-step inversion; see [Reid & Thurston (2014)](https://www.reid-geophys.co.uk/wp-content/uploads/2017/11/Reid-Thurston-2014.pdf). Depths are effective depths below the reference observation plane, not verified ground-depths in a draped survey.
3. **The inherited gate is impossible and circular.** It requires 18/24 strict block wins but only 16 blocks contain proxy truth; eight empty blocks always score zero. The incumbent is itself SGMC-derived. We reproduced the exact frozen proxy/score, disclosed both defects, and did not quietly lower the bar. H4 also loses to H33 on that same instrument, so the gate defect does not rescue it.
4. **Data placement is completed.** The 418,912,844-byte feature raster, template, labels and correct proxy were restored autonomously and hash checked. Same-named sibling SGMC rasters differ; only the pinned GEMSDOE30 copy matches this frozen research instrument. These are documented public mirrors, not an authenticated DrivenData download.
5. **The owner's best is pruning, not new geometry.** H33-B2 is exactly the 40,199-pixel H27-4 base with 2,545 pixels at ≤200 m from the catalogue removed; no new pixels. Reducing false-positive mass can improve DTI without discovering a fault. The **0.2778** file score remains owner-reported; public account rows do not authenticate the TIFF. [Measured analysis](docs/data/h33-measured-analysis.json) / [Explanation](docs/leaderboard-analysis.html).
6. **The leaderboard target changed.** The official page observed on 2026-10-06 showed **0.3345** first, then 0.3262, 0.3222 and 0.3195. This is a dated observation, not a live scrape. We do not claim H4 can beat those scores. [Observation](docs/data/official-observation-20261006.json).
7. **Artifact/site failures are visible.** H33's 216-byte ZIP contains no TIFF; its actual TIFF is sound. Three legacy predictions use `nodata=0`; their stored zeros must be compared, not masked away. Two historical NaN variants also had millions of missing **inside-footprint** values. Their actual finite values match their separately published all-finite companions exactly, with missing cells represented only by zero; the full companions are audited, not silently imputed. Unreadable/unclassified files still fail closed. [Forensic classification](docs/data/prior-nonprediction-forensics-20261006.json). Old downloading/prepare and uniqueness paths could report readiness or novelty with missing data; these are hardened.

## Ranked new hypotheses, before implementation

Four hypotheses were registered before H4 code/scoring: **H4** offset-aware rank-adaptive contact Euler (medium cost, implemented); **H5** two-edge finite density-step inversion (high); **H6** depth-cloud plane geometry / dip localization (high; actual ground projection blocked pending survey datum); **H7** TMI/RTP representation stability without double-counting correlated channels (medium). Expected improvements are qualitative and uncertain, not invented numerical forecasts.

Read the [unchanged detailed registration](docs/research/h4-preregistration-20261006.md), [ranked hypotheses page](docs/hypotheses.html), [research source/claim ledger](docs/research/scientific-library-20261006.md), and [all 44 user-listed projects / 48 supplied scores](docs/leaderboard.html). Euler theory itself is not new science; the artifact is independently generated and novel within the fully specified public corpus.

## Reproduce on CPU

Python 3.11 was used. No GPU, trained neural network, password or manual data placement is required for this Euler pipeline. It is **not** the previously claimed ready-to-train U-Net pipeline.

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-research.txt
python -m pip install -e .
bash scripts/download_competition_data.sh
python scripts/prepare_data.py
python scripts/fetch_prior_cache.py --inventory docs/data/prior-inventory-20261006.json --cache data/prior
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 python scripts/run_contact_euler.py
OPENBLAS_NUM_THREADS=1 python scripts/audit_contact_candidate.py
python scripts/make_contact_previews.py
python scripts/build_contact_site.py
python -m pytest -q
python scripts/build_site.py --check
```

Generation first writes to ignored `work/`, reads **no prior prediction or holdout raster**, and does not use a gradient threshold. Only the separate audit runner can publish a research file after exact-format and full raw-output novelty checks. It evaluates a fixed method without a holdout parameter search and never uploads to the competition. For a **new** experiment, refresh the corpus with `scripts/refresh_prior_inventory.py`, then `scripts/audit_prior_history.py --refresh-heads`, review the exclusions, and audit every byte before claiming novelty; the reproducibility command deliberately uses the frozen dated corpus.

Large inputs/prior caches stay out of Git. Only the deliverable TIFF, compressed cloud, compact evidence and actual-data previews are committed. See [three-pass review](docs/reports/review-20261006.md) and [reproduction receipt](docs/data/h4-reproduction.json). The scheduled Pages workflow refreshes **USGS/GDR context only**, exposes stale/error states, and never changes model values or submission eligibility. Automated DrivenData monitoring is not enabled under its [Terms of Use](https://www.drivendata.org/termsofuse/).

## What remains blocked / next session

- **No independent new-fault holdout or hidden organizer labels.** Register a feasible, source-independent whole-fault-system validation protocol before another model search. Do not convert the circular 0.8359 proxy score into a leaderboard forecast.
- **No DrivenData account session / enrollment here.** No authenticated official download or competition submission was made. No credentials are requested. Hash-pinned public mirrors solved data placement, not source authentication.
- **Observation-surface uncertainty.** Audit official GeoDAWN drape/altitude metadata and native gravity resolution before a ground-depth or dip-to-surface projection.
- **Fault versus contact versus resource.** Euler clusters can be unfaulted lithologic boundaries. Finite-step modeling and independent structural/thermal evidence are needed. The competition target is faults, not vents or commercial reservoirs.
- **Score improvement is unproved.** H4 is a valid, genuinely new *research submission-format candidate*, not a winning candidate. The next experiment must beat a credible independent comparator before spending a slot.

## Concurrent mainline work was preserved, not overwritten

While this branch was running, PRs #7, #8 and #9 landed on main. They are integrated with their source, tests, receipts and downloadable artifacts intact. The current **contact-offset H4** is not the **H4-A contact-network** arm, **H7 gravity-context** arm, or **H40 dotted-emission** arm. H7 in our unchanged preregistration means the still-unimplemented TMI/RTP stability idea, a numbering collision—not the merged H7 experiment.

- **H7 gravity-context:** [separate negative-result report](docs/reports/validation-h7-20261006.md); its 0.000265 proxy result uses **643cbe…**, not H4’s frozen **26d142…** raster. Do not compare those native results as if the labels were identical.
- **H40:** [original artifact/receipt](docs/downloads/gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-audit.json); its claimed promotion estimates were subsequently withdrawn. Both original twins remain research-only, not this H4 file.
- **H4-A:** [closed negative result](docs/reports/h4a-negative-result-20261006.md); its raw raster was withdrawn into another session’s ignored workspace and is not available to this checkout. Do not invent its pixels or pretend it was raw-correlated.
- The [16-pair instrument analysis](docs/data/instrument-audit-20261006.json) and [attribution retraction](docs/data/unauthenticated_attribution_audit.json) disagree in how they describe score evidence. We retain both, but **do not certify those file/score pairings as organizer-authenticated**, infer a hidden-label count from them, or adopt a rank-correlation calibration as a valid promotion gate.

**No candidate is cleared today.** The inherited instrument is **retired for promotion**, but its frozen computation/result remains reproducible. The former missing-proxy blocker is resolved: the exact 26d142… version was recovered from GEMSDOE30. The H7/H40 643cbe… profile and its existing `docs/data/input_manifest.json` are retained separately; H4’s commands use `scripts/acquire_data.py` / `docs/data/acquisition-20261006.json` and never silently substitute profiles. [Integration notes](docs/reports/integration-20261006.md) / [incoming README archive](docs/archive/README-main-9987e40-20261006.md).

## History and operating record

The previous H1/H2 feasibility stops and H2-B negative result remain intact: [append-only validation](docs/research/validation.md), [previous README](docs/archive/README-h2b-20261005.md), [historical H2-B page](docs/index-h2b-20261005.html). Archived mainline outputs are educational references, never this session's candidate. The current site is rebuilt from current JSON receipts; past scores are not silently overwritten.

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
