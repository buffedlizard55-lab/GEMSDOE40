# Three-pass verification record — H8 lineament candidate (2026-10-06)

Every line below was executed in this session. Numbers are read from the receipts named
beside them; nothing is asserted without a receipt, and nothing claimed about the hidden
label set.

## Pass 1 — implement and verify

| Step | Command | Result |
|---|---|---|
| Solver control on an exact analytic contact | `python work/euler_control.py` | windows 9/15/25 accept 585/912/1,417 solutions; median lateral error 0.6/0.3/1.0 m; depth 500 m recovered exactly; σ = 1 → 4 solutions (72 m median lateral error); σ = 5 → 0 |
| H8 generation | `OPENBLAS_NUM_THREADS=2 python scripts/run_h8_lineament.py` | exit 0; 40,000 dots; continuous values 0.383–1.000; receipt `work/h8lineament/h8-lineament-receipt.json` |
| Format gate | `validate_candidate` inside the audit | valid; single band float32; EPSG:32611; shape 3730 × 3292; transform (100, 0, 243350, 0, −100, 4508550); 5,167,373 finite in [0,1]; 7,111,787 NaN outside |
| Novelty gate | `scripts/audit_h8_lineament.py` | 278 comparable priors + 1 unreadable; max abs Pearson 0.0253; max top-mass Jaccard 0.0154; 0 exact duplicates; novel |
| Unit tests | `python -m pytest -q` | 144 passed |

## Pass 2 — review for bugs, missing requirements, wrong assumptions, edge cases

Defects found by review and fixed, with the evidence that they were real:

1. **Direction convention (90° error).** `direction` was measured from +column while
   `anisotropic_kde` projected offsets from +row, so the "along-lineament" smear ran across
   the lineament. Fixed to a single convention (`theta = arctan2(vec[0], vec[1]) mod π`,
   `along = dc·cosθ + dr·sinθ`); covered by `tests/test_h8_lineament.py::test_anisotropic_kde_smears_along_the_lineament`.
2. **Broadcast failure in the KDE.** A leftover per-point `sigmas` line raised
   `ValueError: operands could not be broadcast together with shapes (9,9) (41,)`; rewritten
   as a 12-direction × 5-coherence kernel bank with an `np.add.at` scatter and mass scaled by
   `sigma_along`, so a strung-out train outranks an isotropic blob
   (`tests/test_h8_lineament.py::test_anisotropic_kernel_smears_along_the_train_not_across_it`
   asserts along-axis mass > 5× the across-axis mass on a synthetic NE train).
3. **Solver-control unit bug.** The control compared metres to cells and printed a spurious
   ~1,980 m error; corrected to metres. The corrected output is the row above, and it changes
   a previously published statement in this repository: the solver is *accurate* on model
   data. The correction was published rather than the older, wrong number.
4. **Audit publication crash.** `document["published"]["hard_twin"]` was assigned before
   `document["published"]` existed (`KeyError: 'published'`); fixed by ordering, and the run
   then published both files and the twin's support-identity check.
5. **Portal rejection root cause re-examined.** The family's one observed rejection
   (`Predicted values must be in range [0, 1]`) came from a GeoTIFF written with the integer
   predictor 2 on float data. This release was regenerated with `predictor=1` (no predictor),
   matching the official `sample_submission.tif`, and the TIFF tags were verified with an
   independent struct-level reader (not only rasterio): `Compression = 8`, `Predictor = 1`,
   `SampleFormat = 3`, `BitsPerSample = 32`.
6. **Containment 1.0 explained, not left dangling.** The earlier H8-v1 audit's unexplained
   `max_containment_topk = 1.0` is caused by six dense prior rasters whose positive support is
   the entire footprint (all six have exactly 5,167,373 positive cells; one is a
   positive-everywhere field with minimum 0.00046 and 4,729,578 distinct values).
   `work/h8_v1_containment.json` records all six. Containment is
   therefore reported as a diagnostic only.
7. **Promotion instrument eliminated out-of-sample.** A per-block credit field fitted to 15 of
   the 16 recorded scores by non-negative least squares and tested on the held-out one scores
   **negatively** (LOO Spearman −0.56 to −1.00 across basis grids and assumed hidden counts;
   `work/score_anchored_truth.json`). This is the measured reason the release is not promoted.
8. **Missing artifacts flagged.** The 1 m lidar scarp stack used by an earlier hypothesis is
   absent from this workspace, and `work/h8/h8-field.npy` plus `work/bandcache/` are caches.
   The site and README state the absence instead of quoting an unverifiable number.
9. **Site identity drift.** `scripts/build_contact_site.py` was still asserting the retired H4
   identity while the site served H8; it is now a thin, deterministic wrapper that verifies the
   retained H4 bytes by SHA-256 and then runs the single current generator, so CI cannot pass
   while the served page and the asserted identity disagree.
10. **Non-reproducible gzip container.** The solution cloud was written with
    ``gzip.open()``, which stamps the container with the current mtime: the *content* was
    identical on a re-run but the file hash was not. Fixed with
    ``gzip.GzipFile(..., mtime=0)``; the published cloud was rebuilt (new container hash,
    decompressed bytes proved identical to the previous container) and the whole release now
    reproduces byte-for-byte — continuous TIFF, hard twin and cloud — from an independent
    re-run of the merged runner.
11. **Merge with the published session-2/3 line.** This branch predated mainline sessions 2–3.
    The merge keeps their artifacts (H13 crest-binary, H8 trace-locked depth-KDE, H8-ASA
    analytic-signal), receipts, registers, reports and hashes byte-identical, freezes
    ``src/gemsdoe40/raster.py`` at their recorded source hash (this release writes its
    predictor-1 GeoTIFF with a local writer instead), preserves their manifest as
    ``docs/data/session2-artifacts-20261006.json``, retires their page generator to
    ``scripts/retired/build_h13_site_session2.py``, and renumbers this session's new proposals
    to **H14–H17** because H9–H11 and H12–H13 were already used and closed by them.

## Pass 3 — re-check against the original request

| Requirement | Where it is satisfied |
|---|---|
| Unique single-band float32 GeoTIFF, EPSG:32611, 100 m, exact sample shape/transform, NaN only outside, values in [0,1] | audit §Format; `tests/test_current_release.py` re-reads the published bytes |
| Euler deconvolution depth-clustering (Reid et al. 1990) over magnetic + gravity, fault-like contact SI, depth-labelled solution cloud — not a gradient threshold or edge map | `scripts/run_h8_lineament.py`, `src/gemsdoe40/h8_*.py`; cloud published as `docs/downloads/h8-lineament-solutions.csv.gz` |
| KDE of solution density per pixel, weighted so tight shallow clusters outrank scattered/deep ones, normalised to [0,1] | weighting chain in the generation receipt; `value_fraction_at_one = 0.5`, `distinct_values = 19,923` |
| Hash and correlate against every prior submission before download; refuse to call it new if near-duplicate | audit novelty gate (refuses on exact duplicate or threshold breach); coverage caveat printed on the site and in the audit |
| Easy-to-download TIF at the very top of the site + unique name + short note | hero download on `docs/index.html`, root `index.html` alias, `docs/executive-summary.html` with copy buttons; tracking name `GEMSDOE40-H8-LINEAMENT-785c4f5d5ce1`; note 172 characters |
| Executive-summary subpage explaining exactly how to make a submission | `docs/executive-summary.html`, five numbered steps with the official URLs, format contract and predictor explanation |
| GitHub Pages site, clean and organised, official verified sources | `docs/` + `.github/workflows/pages.yml`; `docs/sources.html` with competition, staff-clarification and geophysics sources |
| 3–5 ranked new hypotheses with layer, physical signature, why-missing, difference, cost, validation | `docs/hypotheses.html` and `docs/research/h8-preregistration-20261006.md`: H14, H8, H15, H16, H17 |
| Full prompt retained in the README and re-read each session | README retains the complete brief verbatim between its markers; `scripts/build_site.py --check` enforces the match |
| No scheduled scraping, no passwords, no automatic submission | feed scope is official USGS/GDR context only with freshness/errors published; no credentials anywhere |
| Three passes; PR and merge; next-session work and blockers | this record; PR body; README "Limitations that remain in the way" |

## Standing negative results (do not re-litigate without new evidence)

* Every local proxy fails to rank the 16 recorded scores (ρ ≤ 0.68); none reaches the ρ ≥ 0.8
  the promotion rule requires, so no local number can justify a weekly slot.
* The counterfactual "off-catalogue SGMC faults are found at catalogue tip extensions" is
  **not** supported: the extension corridor carries 0.0162 credit per off-catalogue pixel,
  0.22× a uniform-random footprint mask of identical size (`work/h9_extension_probe.json`).
  H14 therefore requires the geophysical strike test, not geometry alone.
* The family's best score is pruning, not discovery; adding mass at catalogue-adjacent
  distances is expensive under α = 0.2 / β = 0.8.
