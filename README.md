# GEMSDOE40 — fault-prediction research for the DOE GEMS Prize

> **Current decision: HOLD — DO NOT SUBMIT.** The H2-B research GeoTIFF is unique and passes the local sample-grid/range audit, but its owner-derived SGMC proxy result fails the frozen promotion gate by a wide margin. No weekly submission slot has been used, no DrivenData upload has been made, and no organizer score is known.

**Official target:** geological faults that may indicate geothermal resources—not geothermal vents, hot springs, or reservoir locations. The challenge's private labels and leaderboard are not accessible in this session. This repository is an auditable research and validation system, not a claim of a winning prediction.

## One-click research artifact

- [Download the H2-B research GeoTIFF](docs/downloads/gemsdoe40-h2b-tmi-euler-natural-support-20261005.tif) — **research-only; do not submit**.
- Grid/format: one-band `float32`, EPSG:32611, 100 m, 3,730 × 3,292, exact sample transform/footprint, finite `[0,1]` in-footprint values and NaN outside.
- 17,425 positive cells (0.337% of the 5,167,373-cell footprint); SHA-256: `02486eaa491da2d8ebc6b1cd53ed2bdc8e995f316be1110a64eec484779b5a68`.
- It differs from the 268 same-grid prior rasters under every frozen raw/canonical-hash, correlation, support-Jaccard, and top-budget near-duplicate check. **Uniqueness is not evidence of quality.**
- Suggested DrivenData name and note are intentionally **not provided as upload instructions** while the result is on hold. The site offers a clearly labeled research-only metadata string for record-keeping, not submission.

The full result is on the [project site](docs/index.html), with the [executive summary and future upload guide](docs/executive-summary.html), the [hypothesis register](docs/research/hypotheses.md), and the [H2-B spatial-proxy report](docs/reports/validation-h2b-20261005.md).

## Current evidence (2026-10-05 UTC)

| Stage | Result | Interpretation |
|---|---:|---|
| H1: joint magnetic–gravity Euler depth-consensus KDE | 2,305 accepted magnetic solutions; **0 gravity** solutions; 0 pairs | Preregistered feasibility stop. No candidate or holdout score. H1 was not relaxed or backfilled. |
| H2: four-height TMI Euler persistence, fixed 45,962-cell budget | 17,425 positive KDE cells | Preregistered fixed-budget failure; it aborted rather than inventing support. |
| H2-B: support-limited H2 research raster | 17,425 cells; format pass; uniqueness pass | Append-only support rule registered before proxy scoring. Its proxy result fails promotion. |
| H2-B proxy DTI | **0.006313** | Owner-derived SGMC mirror only; not hidden truth or an official score. |
| H3 ComCat source audit | 34,716 events advertise focal mechanisms; 238 at magnitude ≥4; one detail verified | Usable-plane count and spread across 8+ blocks unverified; H3 remains blocked and unimplemented. |
| Frozen all-prior incumbent | **0.835907** | Same proxy instrument; highest format-eligible same-grid prior raster. The artifact is explicitly SGMC-gap-derived, so this is a circular proxy ceiling—not an independent estimate of challenge generalization. The frozen gate nevertheless retains the all-prior comparator. |
| Promotion | **0/24 blocks won; HOLD** | No weekly slot used; do not submit. |
| Live official leaderboard snapshot | #1 0.3262; #2 0.3222; #3 0.3220; #4 0.3195 | Read from the official board on 2026-10-05. These account-level scores do not identify a TIFF unless the organizer links one. |

The SGMC-derived raster is an owner-derived mirror whose provenance is not independently verified here. It is only a proxy instrument. The best proxy-scoring prior is called `gapfinder-v2-sgmc-gap`; the originating project documents that this candidate used SGMC-gap pixels. Its high proxy DTI is therefore circular. It remains in the frozen all-prior comparison, and the result is disclosed rather than silently removed or relabeled. See [validation report](docs/reports/validation-h2b-20261005.md).

## Project charter — the brief retained for future work

The original project request, captured here as the durable working brief, is to build a unique, scientifically grounded DOE GEMS fault-prediction GeoTIFF together with an auditable research, validation, and project-site system. It calls for:

1. Predict **geological faults that may indicate geothermal resources**, not vents or reservoirs; pursue leaderboard performance without treating any proxy as hidden truth.
2. Preregister three to five ranked hypotheses, implement the scientifically plausible leader, test it on spatially blocked holdouts before using any of the three weekly submission slots, and never invent a score or submit an unvalidated hypothesis.
3. Follow Euler depth-deconvolution/cloud-KDE methods if their physical inputs and feasibility gates support them. Do not replace Euler with a gradient threshold and describe it as Euler. Preserve negative outcomes; do not loosen a frozen gate after seeing its result.
4. Generate a genuinely distinct TIFF, compare raw and canonical hashes, correlations, and support overlap against prior outputs, and disclose any similarity. Copy prior outputs only for education, never as a “new” submission.
5. Match the organizer's sample grid and required GeoTIFF semantics exactly; explicitly prevent the portal's “Predicted values must be in range [0, 1]” error by auditing finite in-footprint values and outside-footprint nodata.
6. Use verified official sources where possible, link evidence, flag source irregularities and access limits, and do not claim official scores, file attribution, or a PR/merge that has not happened.
7. Build a user-friendly GitHub Pages project with a prominent executive summary, one-click GeoTIFF download, current dated leaderboard feed, submission guidance/name/comment when a candidate is genuinely eligible, plus this brief in the README.
8. Keep the operating principles visible: *Maximize P(Win)* and *Own the Outcome*.

**Access boundary:** official challenge data download/submission requires an authenticated account and enrollment. The available feature/sample/label rasters came through a public sibling-repository bridge; the SGMC proxy is an owner-derived mirror. No private test labels or organizer score are present. No manual input has been supplied in this session.

## Research plan and validation protocol

The append-only register is [`docs/research/hypotheses.md`](docs/research/hypotheses.md). It contains three ranked hypotheses, H1's fixed implementation/gates, an H2 lock appended after H1 stopped, and the separate H2-B natural-support rule. The H2-B amendment was entered before its proxy holdout score was read. It did not rewrite or rescue H1/H2.

The frozen promotion gate is intentionally strict: on the same four-by-six spatial blocks and three-cell guards, strictly beat the best format-eligible, comparable same-grid prior raster; win at least 18 of 24 blocks; and lose no block by more than 0.005. Also pass exact-grid/range/nodata validation and the prior-raster near-duplicate audit. A failed gate means **HOLD; DO NOT SUBMIT**. All-prior proxy scores remain visible even when a prior output has circular SGMC provenance; that comparator limitation is explicitly reported.

### Reproduce (after obtaining the required inputs)

Python 3.10+; install with `python -m pip install -e '.[test]'`. The large TIFF inputs and prior-raster cache are external artifacts and are not committed. This run used NumPy 2.4.6, SciPy 1.17.1, Rasterio 1.4.4. Set the paths below to local copies of the pinned inputs and to the verified 270-blob prior cache:

```bash
python -m pytest -q
python scripts/fetch_prior_cache.py --inventory docs/data/prior_raster_inventory.json --cache /path/to/prior-cache
PYTHONPATH=src python scripts/run_research.py \
  --hypothesis H2B \
  --features /path/to/training_features.tif \
  --sample /path/to/example_submission.tif \
  --labels /path/to/existing_faults.tif \
  --proxy /path/to/owner-derived-sgmc-proxy.tif \
  --prior-cache /path/to/prior-cache \
  --inventory docs/data/prior_raster_inventory.json
```

The final command is CPU-heavy (about 7.6 minutes in this sandbox), writes the candidate under `docs/downloads/`, and regenerates validation/format/uniqueness JSON. It **does not submit** anything. H1 and strict H2 have feasibility stops; H2-B is the only completed proxy-scored candidate. The proxy SHA-256 is pinned in the reports and input manifest.

## Official and scientific references

- [DrivenData challenge, task, metric, and submission instructions](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/)
- [Official public leaderboard](https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/) — dated snapshot in [`docs/data/feed.json`](docs/data/feed.json).
- [DOE/NLR September 2026 competition rules PDF](https://docs.nlr.gov/docs/fy26osti/96647.pdf)
- [DrivenData staff clarification: exact-pixel known-fault mask](https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516/2)
- [USGS State Geologic Map Compilation metadata](https://mrdata.usgs.gov/geology/state/USGS_SGMC_Metadata.xml) and [Data Series 1052](https://pubs.usgs.gov/ds/1052/): official source context for the vector compilation; this does **not** validate the owner-derived raster mirror used here.
- Reid & Thurston (2014), [structural-index guidance for gravity and magnetic interpretation](https://www.reid-geophys.co.uk/wp-content/uploads/2017/11/Reid-Thurston-2014.pdf).
- [USGS ComCat FDSN Event Web Service](https://earthquake.usgs.gov/fdsnws/event/1/) — official H3 source; preliminary product counts and one detail example are in [`h3_source_audit.json`](docs/data/h3_source_audit.json), but the registered usable-plane/spatial-coverage gate is still unmet.
- [GEMSDOE28 Euler research](https://buffedlizard55-lab.github.io/GEMSDOE28/) — prior art, not an official source; H1/H2 novelty claims are narrow and explicitly bounded.

## Repository map

- `docs/index.html` — project status and one-click research artifact.
- `docs/executive-summary.html` — decision, audit trail, and honest future submission instructions.
- `docs/research/hypotheses.md` — append-only preregistration and implementation locks.
- `docs/research/source-register.md` — verified links, provenance limits, and data irregularities.
- `docs/reports/` — human-readable feasibility and validation reports.
- `docs/data/` — pinned input manifest, prior inventory and full metadata audit, exact receipts, and machine-readable feed/results.
- `src/gemsdoe40/` — metric, Euler generators, holdout scoring, raster validation, uniqueness audit.
- `scripts/` — reproducible cache and research runners.
- `tests/` — numerical, format, and persistence regression tests.

**Status discipline:** no file here has an organizer score. No submission has been made. A public account-level leaderboard row does not prove which TIFF earned it. Any future PR, merge, page deployment, or submission must be verified from its actual service response before being described as completed.
