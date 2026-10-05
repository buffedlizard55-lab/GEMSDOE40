# GEMSDOE40 — Euler deconvolution depth-clustering submission system

**Competition:** DOE GEMS Prize (DrivenData #306) · [overview](https://www.drivendata.org/competitions/306/competition-doe-gems/) ·
[problem page](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/) ·
[rules PDF](https://docs.nlr.gov/docs/fy26osti/96647.pdf) · **live site:** <https://buffedlizard55-lab.github.io/GEMSDOE40/>

> ### Our Core Values — kept central to every decision in this repository
> **Maximize P(Win).** *"Maximize the Probability of Winning"* is our decision-making framework. In
> every decision we weigh tradeoffs, assess risk, and choose the path that maximizes the probability
> that we win this prize. We set aside our emotions and make tough decisions in order to maximize
> P(Win). It frees us from constraints and clarifies that we must put the outcome first.
>
> **Own the Outcome.** We own results end to end — not just our individual slice of the work. When
> problems arise and we have the means to act, we act without waiting for permission or assignment.
> We treat failure and success as signals and use them to improve. We stay accountable to the final
> outcome.
>
> Applied here: *Maximize P(Win)* is why the session ends with a measured **upload / do-not-upload**
> verdict instead of a file dump, and why a hypothesis that failed (raw Euler centroids, Euler as a
> gate) is published as failed. *Own the Outcome* is why every number on the site is regenerated from
> bytes by a script in `scripts/`, and why the one defect found in our own first attempt (window
> windows in the Fourier taper zone scoring spurious 20 m solutions) is recorded in the test suite
> rather than quietly patched.

**Standing brief for every session:** [`docs/prompt.md`](docs/prompt.md) — read it first.

---

## ⬇ ONE-CLICK SUBMISSION FILES

| file | what it is | measured | verdict |
| --- | --- | --- | --- |
| **[`gems40-euler-augmented-incumbent-20261005-ac1c1927.tif`](docs/downloads/gems40-euler-augmented-incumbent-20261005-ac1c1927.tif)** (885 kB, sha256 `54aa8c6cbee0a1f8…`, 60,710 px) | current site-best emission ∪ **every off-catalogue Euler SI=0 depth-cluster crest** (23,056 new px) | LM-calibrated **0.30845** vs the site best's 0.26792 → **+0.04053, better in 4/4 blocked folds**; 18/18 format checks pass | **UPLOAD THIS** — flagged NEAR-DUPLICATE by the novelty audit *by design* (it contains the previous submission's pixels; it is an augmentation, not a new pattern) |
| **[`gems40-euler-si0-depthcluster-crossfamily-20261005-8933d380.tif`](docs/downloads/gems40-euler-si0-depthcluster-crossfamily-20261005-8933d380.tif)** (635 kB, sha256 `7ee6cac2444949a3…`, 20,000 px) | the brief's mandate: magnetic × gravity Euler SI = 0 depth-clustering, KDE → normalised [0, 1] | novelty verdict **NEW** (max Jaccard 0.0049, max |r| 0.0050 vs 16 prior rasters); LM-calibrated 0.11600 | **RESEARCH ONLY** — a genuinely different spatial pattern, but it does not beat the incumbent standalone, so the standing rule forbids spending a slot on it |
| [`…-zeros.tif`](docs/downloads/gems40-euler-si0-depthcluster-crossfamily-20261005-8933d380-zeros.tif) companions | identical inside the footprint, `0` instead of `NaN` outside | all values finite in [0, 1] everywhere | use if the portal rejects NaN |
| [`…-continuous-….tif`](docs/downloads/gems40-euler-si0-depthcluster-crossfamily-continuous-20261005-f3e0ed7a.tif) | the literal reading of the brief: continuous normalised KDE, no crest reduction | LM-calibrated 0.12487, 3.38 M px | **AUDIT ONLY** — the metric charges α per emitted pixel |

Full receipts (12 portal checks each): [`docs/downloads/manifest.json`](docs/downloads/manifest.json).

**How to submit in 30 seconds:** download the first file → open the
[DrivenData submission page](https://www.drivendata.org/competitions/306/competition-doe-gems/submissions/)
→ upload it → paste the note from [`docs/executive-summary.html`](docs/executive-summary.html#3-step-by-step)
→ submit. That page also explains the **"Predicted values must be in range [0, 1]"** error you saw
and how every file here avoids it (finite, in [0, 1] inside the 5,167,373-pixel footprint; NaN only
outside it, exactly like the organizers' own sample submission).

---

## What was done in this session

1. **Read the whole brief and the programme's history.** The score record, the metric algebra and the
   sibling repositories' instruments were read from their published knowledge files
   (`GEMSDOE28/knowledge/43_hypotheses_heatflow_euler.md`,
   `GEMSDOE32/knowledge/02_the_ceiling_and_the_instrument.md`, `GEMSDOE32/src/gems32/live_mirror.py`).
2. **Implemented Euler deconvolution from scratch and verified it against analytic sources.**
   `src/gems40/euler.py` (Fourier derivatives → box-filtered sliding-window least squares → depth
   estimate per window), `tests/test_euler.py`: a monopole (N = 2) and a 2-D contact (N = 0,
   `T = atan(x/z0)`) recover known depths (688 m recovered vs 700 m true for the contact), and the
   vertical-derivative sign convention is asserted by flipping the sign and requiring negative
   depths. No external library was used for the inversion.
3. **Built the depth-aware weighting and the KDE raster** (`cluster_weights`, `kde_field`), ran it on
   all four potential-field layers of the pinned competition stack — `rtp`, `tmi`, `mag_anom`
   (magnetic) and `iso_grav_anom` (gravity) — with SI = 0 (fault/contact, Reid et al. 1990 Table 1).
4. **Verified the decision instrument before using it.** `scripts/validate_candidates.py` re-derives
   the LM instrument (blocked quadrants, off-catalogue SGMC truth, prevalence-calibrated DTI) and
   reproduces the published live orderings that matter at this operating point:
   0.2778 > 0.2708 > 0.2600 > 0.2477.
5. **Measured the hypotheses and published the failures.** Raw Euler KDE mass adds credit at
   0.040–0.053 credit/px — *below* the metric's break-even τ = 0.0588. Crest-thinning the same cloud
   raises it to **0.066–0.085 = 1.13–1.44 × τ**, positive in 4/4 folds. Euler as a *gate* on the
   incumbent is **refuted** (crest support and incumbent support are nearly disjoint, Jaccard ≈ 0.005).
6. **Produced the unique artifact the brief demands and proved it is new**:
   `scripts/audit_artifacts.py` hashes and correlates every shipped raster against all 16–17 prior
   rasters on disk under a rule fixed *before* the numbers were seen (near-duplicate ⇔ Jaccard > 0.50
   or |Pearson| > 0.80). The pure Euler artifact: max Jaccard **0.0049**, max |r| **0.0050** →
   verdict **NEW**. The augmented artifact is honestly reported as NEAR-DUPLICATE because it embeds
   the previous submission.
7. **Published the site** (`docs/`, GitHub Pages) with the executive summary, method, ranked
   hypotheses, leaderboard analysis and evidence tables — all generated by `scripts/build_site.py`
   from the JSON evidence, so nothing on it can drift from the measurements.

## Honest verdict and expected value

* The unique Euler artifact is **new information** (its support is essentially disjoint from every
  prior submission) — this is the first genuinely independent channel this programme has added since
  the LiDAR-scarp ridge family.
* Its mass is worth **1.13–1.44 ×** what the metric charges for it, so adding it to the incumbent is
  expected to **raise** the score. Conservative projection (live-calibrated |G| = 12,226 and the
  measured credit/px): **≈ 0.288**; instrument-based projection (+0.0405): **≈ 0.318**. The LM
  instrument is documented to over-reward denser emissions, so plan with the conservative figure.
  **Neither figure reaches 0.3195**, and this repository says so plainly.
* The remaining gap to #1 is a *detector-quality* gap, not an emission-rule gap. The two cheapest
  next tests (multi-window Euler conjunction; depth-banded emission) are named on the
  [hypotheses page](docs/hypotheses.html) with their cost and their data needs.

## Reproduce

```bash
bash scripts/fetch_data.sh            # sha256-verified mirrors of the competition rasters
python scripts/inspect_data.py        # re-derive every data claim from the bytes
python scripts/run_euler.py --layer rtp          # also: tmi, mag_anom, iso_grav_anom
python scripts/validate_candidates.py             # instrument verification + candidate sweep
python scripts/union_test.py                      # marginal value on the incumbent (break-even test)
python scripts/nms_refine_test.py                 # crest refinement
python scripts/euler_gate_test.py                 # Euler as a prune/gate (refuted)
python scripts/build_submission.py                # write the artifacts + receipts
python scripts/audit_artifacts.py                 # validation + hash/correlation novelty audit
python scripts/build_site.py                      # regenerate docs/
python -m pytest tests -q                         # 10 tests, incl. analytic Euler verification
```

Requirements: `pip install -r requirements.txt` (numpy, scipy, rasterio, pytest). CPU only; the whole
pipeline runs in ≈2 minutes on 2 cores.

## Repository map

| path | purpose |
| --- | --- |
| `docs/prompt.md` | the standing brief (read every session) |
| `src/gems40/euler.py` | Fourier derivatives, sliding-window Euler inversion, depth-aware clustering, KDE |
| `src/gems40/metric.py` | the official distance-weighted Tversky index (soft and binary) + break-even identity |
| `src/gems40/instrument.py` | the blocked holdout instruments (LM off-catalogue, catalogue-component) |
| `src/gems40/grid.py`, `layers.py`, `pins.py` | submission I/O + the 12 portal checks, official band inventory, sha256 pins |
| `scripts/` | fetch → inspect → invert → validate → build → audit → site |
| `data/evidence/` | every measurement as JSON (validation, novelty, sweeps, gate test, format) |
| `docs/` | the GitHub Pages site, generated from `data/evidence/` |

---

## Other artifacts already in this repository (parallel track on `main`)

An earlier scaffold (merged on `main` before this branch) produced a further Euler SI = 0 depth-KDE
artifact, kept here rather than deleted so both tracks stay auditable:

| file | size | that track's own reported measurement |
| --- | --- | --- |
| [`docs/downloads/gemsdoe40-euler-si0-depthkde-20261005T230506Z-5bc279cd-zeros.tif`](docs/downloads/gemsdoe40-euler-si0-depthkde-20261005T230506Z-5bc279cd-zeros.tif) (+ `-nan`, `.zip`, `GEMSDOE40-submission.tif/.zip`) | 340 kB, 45,784 px | uniqueness vs 9 priors: worst Pearson 0.0072, worst Jaccard 0.0054 → new; but on the **catalogue** holdout it reports mean blocked DTI 0.0150, i.e. −0.0674 vs a uniform-random control and −0.0066 vs a gradient baseline (`evidence/last_run.json`) |

That catalogue holdout is the instrument this repository deliberately does **not** decide with: its
truth *is* the published catalogue, and the organisers mask catalogue pixels out of scoring, so it
rewards exactly the mass the live scorer throws away (measured in the sibling repos as ρ = +0.14
against live scores, versus ρ = +0.53 for the off-catalogue LM instrument used here). Both tracks
agree on the science (Euler SI = 0 depth-KDE, catalogue-buffered, `[0, 1]`, NaN outside the
footprint); they differ in what they measured, and both measurements are quoted above rather than
harmonised.

## Limitations (stated, not hidden)

* No DrivenData login and no access to `drivendata.org`, `dropbox.com`, `gdr.openei.org`,
  `docs.nlr.gov` or any `*.github.io` page from this sandbox (all HTTP 000). Competition data and
  prior artifacts were therefore taken from **byte-identical GitHub mirrors** and verified against
  sha256 pins; every live score quoted is **owner-reported**, not an organizer receipt.
* The LM instrument is a proxy: it reproduces 4/6 published live orderings (including the three that
  matter at this operating point) and over-rewards mass beyond ≈120 k emitted pixels. The recommended
  artifact emits 60,710 px, inside that domain.
* The Euler inversion uses one window scale per run; multi-scale conjunction (H40-4) is specified but
  not yet measured.
