# GEMSDOE40 — fault-prediction research for the DOE GEMS Prize

## Permanent owner brief — read this before every session

**Mission:** Maximize P(Win) and Own the Outcome while building an auditable, useful system for predicting *geological faults relevant to geothermal resources*. The target is faults—not geothermal vents or reservoir locations.

### Required scientific and delivery work

1. Generate an original, non-copied, single-band GeoTIFF. For Euler hypotheses, solve potential-field source locations and depths and project a **weighted source-solution cloud KDE**; do not substitute gradient thresholding and call it Euler. Explain structural-index and geologic confounds.
2. Before implementing a new strategy, preregister **3–5 distinct, ranked geological hypotheses**. For every one specify exact layers, a physical signature, why it may find faults absent from USGS/INGENIOUS, how it differs from repository methods, expected DTI direction/rank, implementation cost, and external-data/access needs. Do not call an externally dependent idea viable until free official access is checked.
3. Validate the leader on spatially blocked holdout data before spending any of the official weekly feedback slots. No slot unless it beats the then-current comparable holdout best under frozen gates. A proxy is not organizer truth; never label a proxy result an official score.
4. Compare raw hashes and canonical pixel values against every inventoried prior output. Keep source refs, Git blobs, hashes, correlation and support-overlap measures. Do not copy prior submissions as a new candidate.
5. Re-open the output TIFF and verify the actual sample grid: one band, float32, EPSG:32611, 100 m, exact dimensions/transform/footprint, finite `[0,1]` values inside and NaN outside. Prevent portal range/nodata errors, but treat format as a gate—not a quality result.
6. Use official, trusted sources; link them for manual review; flag source irregularities and access limits; verify claims; do not hallucinate scores, provenance, or file attribution.
7. Make the TIFF download obvious on the project site. Maintain an executive-summary guide that names the exact upload steps, but give a submission name/comment only after a candidate is cleared. Research-only files must say **DO NOT SUBMIT**.
8. Run tests and independent file audits. If a PR is requested, open it from this session's working branch; report a merge only after verifying GitHub's actual response. Never request credentials in chat.

**Verified starting links:** [DrivenData task, metric and format](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/); [DOE/NLR 2026 rules](https://docs.nlr.gov/docs/fy26osti/96647.pdf); [official public leaderboard](https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/); [Reid et al. (1990) Euler method](https://doi.org/10.1190/1.1442774); [Reid & Thurston (2014) structural-index caveats](https://www.reid-geophys.co.uk/wp-content/uploads/2017/11/Reid-Thurston-2014.pdf); [structured standing brief](docs/prompt.md); [2026-10-06 preregistration](docs/research/preregistered-hypotheses-20261006.md).

## Current decision — 2026-10-06 UTC

> **HOLD — DO NOT SUBMIT H7.** H7 passes the exact sample-grid/format audit and is distinct from all 304 comparable public prior TIFF blobs, but scores only **0.0002650893** on the current owner-derived SGMC proxy and wins **0/24** guarded spatial blocks against the frozen all-prior incumbent. No weekly slot was used. No DrivenData upload or organizer score is known.

### Current research TIFF

- [Download H7 research-only GeoTIFF](docs/downloads/gemsdoe40-h7-rtp-euler-gravity-context-3d-kde-20261006-998f660f.tif) — **not cleared; do not submit**.
- Method: 3,237 accepted magnetic RTP SI=0 Euler solutions, 430 supported 3-D source voxels, a preregistered `[0.5,1]` soft gravity-horizontal-gradient context factor, and weighted source-cloud KDE. It is not an edge threshold and does not claim raw-gravity Euler solutions.
- Output: one float32 band, EPSG:32611, 3,292 × 3,730, 100 m, exact sample transform, finite `[0,1]` in-footprint values, NaN outside; 28,080 positive in-footprint pixels.
- Raw SHA-256: `e229aa9af26018bc80b0cca880d484a5b8362c26ad79ca362423953c32e66597`.
- Canonical pixel SHA-256: `998f660fcbf3b1a31600bd1e5e56a4d3b4c32d57c8fd8af36b3af793226a7d7c`.
- Uniqueness: PASS against 304 exact-grid prior outputs (0 near duplicates, fixed 45,962-cell top-support budget). Public-corpus uniqueness is not a performance result or a guarantee against private/unseen artifacts.
- Holdout: current proxy DTI `0.0002650893`; highest all-prior comparator `0.8257222257` (SGMC-gap-derived and circular on this proxy); pooled delta `−0.8254571364`; 0/24 block wins; worst block delta `−0.9421647993`. The highest non-circular prior scored `0.7618701028`. These are **local diagnostics only**, not organizer scores.
- Full report: [`docs/reports/validation-h7-20261006.md`](docs/reports/validation-h7-20261006.md). The executive summary provides a future-only upload checklist; there is no current upload name or note because no TIFF is cleared.

### Frozen H7 receipts

- Format: [`docs/data/format-receipt-h7-20261006.json`](docs/data/format-receipt-h7-20261006.json).
- Construction without proxy read: [`docs/data/h7-build-20261006.json`](docs/data/h7-build-20261006.json).
- 24-block result and conservative promotion gate: [`docs/data/validation-h7-20261006.json`](docs/data/validation-h7-20261006.json).
- Scores for all 304 exact-grid priors: [`docs/data/prior-holdout-scores-h7-20261006.json`](docs/data/prior-holdout-scores-h7-20261006.json).
- Raw/canonical pixel audit: [`docs/data/uniqueness-audit-h7-20261006.json`](docs/data/uniqueness-audit-h7-20261006.json).
- Complete pinned prior inventory: [`docs/data/prior_raster_inventory-20261006.json`](docs/data/prior_raster_inventory-20261006.json).

## Candidate and evidence ledger

| Work | Outcome | Decision and caveat |
|---|---|---|
| H1 / earlier joint potential-field Euler | Feasibility stop | No accepted gravity cloud; no candidate. Do not relax the gates. |
| H2 / fixed-budget TMI persistence | Feasibility stop | Only 17,425 positive cells against the fixed 45,962-cell budget; it aborted rather than fabricate/backfill support. |
| H2-B / natural-support TMI Euler KDE | Historical research candidate | Format and uniqueness passed, but old-proxy DTI `0.006313285` failed the gate at 0/24. Same candidate on the **current** proxy scored `0.006223113`; these results are tied to separate hashes. Still **DO NOT SUBMIT**. See [`validation-h2b-20261005.md`](docs/reports/validation-h2b-20261005.md). |
| H3 / ComCat focal-mechanism constraint | Blocked | Need usable nodal-plane coverage and spatial spread across holdout blocks. USGS ComCat is official/free, but the registered coverage gate is unverified; no H3 implementation. |
| H4 / raw gravity SI=−1 + magnetic Euler depth consensus | Feasibility stop before proxy | 3,237 accepted magnetic sources; 316,339 raw-gravity windows screened, 24 conditioned, **zero** within locked 100–2,500 m depth gate. No candidate; gates were not relaxed. |
| H5 / conductivity front × basement-depth break | Registered, untried | Existing `cond_surf`, `depth_to_base_surf`, RTP and raw gravity layers; no new external data. Geological and alteration confounds remain. |
| H6 / geodetic strain-release relay | Registered, untried | Existing strain and seismic-context layers; no new external data. Spatial validation required before any promotion. |
| H7 / gravity-context magnetic Euler 3-D KDE | Format + uniqueness pass; holdout failure | DTI `0.0002650893`, 0/24 block wins. **HOLD — RESEARCH ONLY — DO NOT SUBMIT.** |

The 2026-10-06 registration contains four ranked, distinct H4/H5/H6/H7 hypotheses; H7 was registered as a fallback before H7 code and output, after H4's locked feasibility stop. Full physical signatures, layer names, novelty boundaries, expected improvement/cost, and locked tests are in [`docs/research/preregistered-hypotheses-20261006.md`](docs/research/preregistered-hypotheses-20261006.md). The earlier H1/H2/H3 records remain append-only at [`docs/research/hypotheses.md`](docs/research/hypotheses.md).

## Proxy, leaderboard, and provenance guardrails

The proxy discrepancy was resolved as two distinct, pinned owner-mirror versions. The current file is SHA-256 `643cbe992ef4ba37588fb469163ed8291e3ceb23d6c1f78a3cfaa462430c2da0`; the historical H2-B file is `26d142c4c93282cd94f6950ab96f22aeff59fbbea523d43d662e76fa1b161b5c`. They differ in 1,450 cells (Pearson `0.991230`; positive-support Jaccard `0.982655`) and must not be merged or presented as organizer ground truth. See the [dated pin reconciliation](docs/reports/proxy-pin-reconciliation-20261006.md), [current input manifest](docs/data/input_manifest.json), and [preserved legacy manifest](docs/data/input_manifest-20261005-h2b-legacy.json).

The official public leaderboard snapshot retrieved 2026-10-06 lists #1 at `0.3345`; `0.3195` is rank 4. The account-level `extradr19` row is rank 13 at `0.2778` (ten submissions in the fetched snapshot). GEMSDOE32 labels H33-2-B2's `0.2747` as **projected / UNSCORED**. No file-level receipt establishes that the account row belongs to H33-2-B2; keep the claims separate. Snapshot: [`docs/data/feed-20261006.json`](docs/data/feed-20261006.json); explanation: [`docs/leaderboard-analysis.html`](docs/leaderboard-analysis.html).

## Prior-raster inventory refresh

- 55 public owner GEMSDOE repository heads; 55 recursive trees read without truncation.
- 501 current-head TIFF paths explicitly classified; zero unclassified paths.
- 48 submission-looking paths were reviewed row by row and classified as prior submission/prediction outputs based on `submission.tif` / `submission_conformant.tif` names and repository context.
- 306 unique Git blobs across 541 pinned repo/ref/path locations; 27 new unique output blobs added over the prior inventory.
- All 306/306 TIFF blobs fetched and Git blob SHA-1/SHA-256 verified in `/tmp/gemsdoe40-prior-cache` (external temporary cache, not committed).
- H7 holdout comparison: 304 exact sample-grid, single-band priors; 183 also pass the strict local format filter; two cached TIFFs were not exact-grid/single-band comparators.

Read the [inventory refresh report](docs/reports/prior-raster-refresh-20261006.md) and full [pinned inventory](docs/data/prior_raster_inventory-20261006.json). Public repository paths do not prove official upload, score, selection, or account attribution; private/unpublished outputs remain outside the census.

## Input facts and irregularities

`python scripts/inspect_data.py` re-verified the locally available input hashes on 2026-10-06. The sample grid is 3,730 × 3,292, EPSG:32611, 100 m, transform `[100, 0, 243350, 0, -100, 4508550]`, with 5,167,373 footprint pixels and 7,111,787 NaNs outside. The sample contains 60,988 exact known-label positive pixels despite sibling prose describing an all-zero sample; it is used only for geometry/footprint, not truth. See [`data/evidence/submission_format.json`](data/evidence/submission_format.json) and [`docs/data/input_manifest.json`](docs/data/input_manifest.json).

The official competition data/submission pages are login-gated in this environment. Inputs were read from public sibling-repository mirrors, byte-pinned, and labelled with their authentication limitations. Neither successful SHA verification nor official USGS metadata authenticates the proxy TIFF as an organizer-provided label raster. There are no private test labels or organizer receipts here.

## Project site and future upload guide

The [project home](docs/index.html) makes the H7 research TIFF download visible and labels it **DO NOT SUBMIT**. The [executive summary](docs/executive-summary.html) records why no file is cleared and gives a checklist only for a future candidate that passes every gate. No submission name/comment is supplied for H7. The [sources and limitations page](docs/sources.html) links primary sources for manual review and records provenance caveats.

A public account leaderboard score, a model projection, a local proxy score, or a format PASS must never be rewritten as an organizer score. No upload has been made.

## Reproducibility

Python 3.10+; install with `python -m pip install -e '.[test]'`. Large TIFF inputs and the prior cache are external artifacts. In this sandbox, use `.venv/bin/python` because system Python may not have Rasterio.

```bash
.venv/bin/python scripts/inspect_data.py
.venv/bin/pytest -q
.venv/bin/python scripts/validate_h7_candidate.py \
  --candidate docs/downloads/gemsdoe40-h7-rtp-euler-gravity-context-3d-kde-20261006-998f660f.tif \
  --inventory docs/data/prior_raster_inventory-20261006.json \
  --cache /path/to/verified-306-blob-cache
.venv/bin/python scripts/build_site.py --check
```

`validate_h7_candidate.py` enforces explicit current-proxy SHA and verifies feature, label, and sample pins; it checks all prior scores, the strict exact-grid format receipt, the complete-cache requirement, and the fixed 45,962-cell uniqueness budget. It writes audit JSON only and does **not** submit. Do not rerun/tune H7 as a search against the same proxy after observing its failed result.

Tests at this update: **42 passed**. The full H7 validation scan scored 304 same-grid rasters and compared all 306 unique blobs; see the linked dated receipts.

## Remaining work and limitations

1. H7 is not competitive on the only available blocked proxy. Do not spend a weekly slot, and do not alter its gates after seeing the result.
2. H5/H6 remain registered but untried. Any new implementation must be selected without proxy tuning, preserve the original registration, and go through a fresh spatially blocked holdout and full public-prior comparison.
3. H3 needs official ComCat nodal-plane coverage/spatial-spread verification before implementation. External-data ideas require official-free-access checks before being called viable.
4. The local owner-derived proxy is not independently authenticated as organizer truth. Best-prior proxy performance is circular in part, and no result predicts the private test or final expert-expanded round.
5. The prior-raster corpus covers public repository default-branch heads at pinned commits, not private repositories, local submissions, or future outputs.
6. No candidate is cleared, no official score exists for these local files, and no weekly slot or upload has been used.

**Official/source links:** [problem description](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/); [rules PDF](https://docs.nlr.gov/docs/fy26osti/96647.pdf); [public leaderboard](https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/); [USGS SGMC metadata](https://mrdata.usgs.gov/geology/state/USGS_SGMC_Metadata.xml) (source context only); [INGENIOUS regional data](https://gdr.openei.org/submissions/1391); [USGS ComCat API](https://earthquake.usgs.gov/fdsnws/event/1/); [Reid et al. (1990)](https://doi.org/10.1190/1.1442774); [Reid & Thurston (2014)](https://www.reid-geophys.co.uk/wp-content/uploads/2017/11/Reid-Thurston-2014.pdf).
