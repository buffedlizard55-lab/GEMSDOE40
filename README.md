# GEMSDOE40 — Euler deconvolution depth-clustering for the DOE GEMS Prize

**Read this file, then `docs/SESSION_BRIEF.md`, at the start of every session.**

The goal is first place on the DrivenData GEMS Prize leaderboard
(https://www.drivendata.org/competitions/306/competition-doe-gems/).
The site’s job is to make submitting a legal GeoTIFF a one-click download,
and to keep an auditable, up-to-date feed so nobody has to check the
format by hand.

## ⬇ ONE-CLICK SUBMISSION FILE

- **Download:** [docs/downloads/GEMSDOE40-submission.tif](docs/downloads/GEMSDOE40-submission.tif)
- **Zip:** [docs/downloads/GEMSDOE40-submission.zip](docs/downloads/GEMSDOE40-submission.zip)
- **How to submit (5 steps):** [docs/executive-summary.html](docs/executive-summary.html)
- **Site:** [docs/index.html](docs/index.html)

**Unique name:** `GEMSDOE40-euler-si0-depthkde`

**Note to paste (106 / 200 characters):**

```
GEMSDOE40 Euler SI=0 depth-KDE | RTP+isograv win8/12, 45784 px, 0 on cat, bin Δvs|∇RTP| +0.0032 | 5bc279cd
```

Primary file is **all-finite float32**, every value in `[0, 1]`, `nodata`
unset, EPSG:32611, 100 m, 3730×3292. This is the fix for the portal error
`Predicted values must be in range [0, 1]`.

sha256 `5bc279cd55e65733df6d532b6b1c9d348b33720ed9196237371d487042fff195`

---

## Standing orders (from the session prompt — do not skip)

**Core values.** Maximize P(Win). Own the Outcome. Weigh tradeoffs, take
the path that maximises the probability that this project places. Own
results end to end.

**Highest urgency.** Generate a *unique* TIF submission. Do not copy a
previous GEMSDOE submission except to learn from it. Hash and correlate
against every prior raw output and refuse to call it new if it is a
near-duplicate.

**Method required this session.** Euler deconvolution depth-clustering,
built from potential-field depth estimation, **not** gradient
thresholding:

1. Run Euler deconvolution (Reid, Allsop, Granser, Millett & Somerton,
   *Geophysics* 1990, https://doi.org/10.1190/1.1442774) across the
   magnetic and gravity layers with the structural index for a
   fault-like contact (SI = 0).
2. Produce a cloud of depth-labelled solution points, not a single edge
   map.
3. Convert that cloud to a continuous raster by kernel-density estimation
   of solution density per pixel, weighted so tight clusters of shallow,
   mutually-consistent solutions score higher than scattered or deep ones
   (a real near-surface fault produces the former; noise produces the
   latter).
4. Normalise to `[0, 1]`. Write a single-band float32 GeoTIFF in
   EPSG:32611, 100 m, matching the sample submission’s exact shape and
   geotransform. NaN only outside the valid footprint — and, because the
   portal has rejected NaN/sentinel files, also ship an all-finite zeros
   twin as the primary download.

**Before implementing.** Generate 3–5 candidate geological hypotheses we
have not tried yet (layer, signature, why it catches a fault *missing*
from USGS/INGENIOUS, how it differs from prior GEMSDOE code). Rank by
expected DTI vs cost. Validate the top candidate on a spatially-blocked
holdout before spending a weekly slot. See
[docs/hypotheses.html](docs/hypotheses.html).

**No hallucinations.** Work line by line from official sources. Provide
links. Flag irregularities. The site must stay current so we do not
manually re-check format.

**Competition.**
https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/

---

## Why 0.2778 scored, and whether we can beat it

GEMSDOE32’s `h33-h33-2-b2` (owner-claimed **0.2778**) is not a new
detector. It is a 2-pixel catalogue-flank prune of a 0.2708 dotted file
which is itself a re-pack of the 0.2600 `dotted-h19-5-d2-8` surface
(44,090 dots at ~300 m spacing). The metric

```
DTI = T / (0.2 (T + S − M) + 0.8 |G|)
```

(α = 0.2, β = 0.8, triangular kernel R = 300 m; official worked example
TP=3, FP=1.89, FN=2 → 0.60, reproduced in `tests/test_metric.py`) is a
**budget**. Mass that is not the best cover of a hidden-truth pixel
costs 0.2. Thinning a thick surface while keeping its geometry raises
credit per pixel. The public leader on 2026-10-05 is **0.3262**
(nchuzhoy), not 0.3195.

Beating 0.3262 requires higher credit density *on the hidden set*. The
hidden set is faults that are **not** in the USGS/INGENIOUS catalogue, so
catalogue-hugging pixels are worth ~0. Copying H33-2-B2 cannot be “new”
and this repository refuses to ship a near-duplicate (measured Pearson
**0.007**, Jaccard **0.005** vs that file).

Euler depth-clustering is a different physical observable: the *(x, y,
z)* of a contact, not an edge of the field. Tight shallow clusters are
the signature of a real near-surface fault; scattered deep solutions are
noise (Reid et al. 1990, figs. 1 and 5). Whether those contacts coincide
with the hidden expert traces is an empirical question that needs a live
slot. **This file is unscored. We do not claim it beats 0.2778.**

What we *did* validate, on the spatially-blocked / same-count instrument
we have:

- Unmasked Euler cores, same count as themselves, vs the public
  catalogue: binary DTI **0.0258**.
- Same-count |∇RTP| gradient baseline: **0.0225**. Δ **+0.0032**.
  Euler is not gradient thresholding.
- 98.4 % of unmasked cores are off-catalogue; the shipped file has **0**
  on-catalogue pixels (2 px / 200 m buffer).

Catalogue-proxy DTI is low because the method is off-catalogue. That is
the *correct* bias for this contest, and it is why we do not use
catalogue DTI as a promotion gate (GEMSDOE32 measured the same
disagreement between catalogue-holdout and live-anchored truth).

---

## Reproduce

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
# data/training_features.tif must be present (sha256 4371c82e…)
# assemble from data/bridge parts, or DrivenData → data/
bash scripts/download_competition_data.sh
.venv/bin/python scripts/prepare_data.py
PYTHONPATH=src .venv/bin/python scripts/run_euler_depth_cluster.py
PYTHONPATH=src .venv/bin/pytest tests -q
```

No GPU. Euler + KDE on the full GeoDAWN grid is ~60 s on CPU.

---

## Limitations (what this sandbox cannot do)

| Limitation | Effect | What we need |
|---|---|---|
| No DrivenData auth | Cannot auto-download a fresh `training_features.tif` from the data tab; cannot press Submit; cannot read a live score for this file | A logged-in machine, or the sha256-pinned data-bridge (used here) |
| Dropbox / raw.githubusercontent.com TLS failed here | IR-40-TLS-01 | GitHub API `Accept: application/vnd.github.raw` worked and was used |
| Holdout truth = public catalogue | Cannot reward a genuinely new fault (IR-40-PROXY-01) | A live slot, or an off-catalogue proxy population (SGMC was used in GEMSDOE31; not rebuilt here) |
| Weekly slot cap | This repository packages a candidate; it does not spend a slot | Human upload of `GEMSDOE40-submission.tif` |

---

## Layout

```
docs/index.html          ← download is the first thing you see
docs/executive-summary.html
docs/hypotheses.html     ← 5 ranked hypotheses; H40-A shipped
docs/sources.html        ← official links, fetched this session
docs/downloads/          ← the GeoTIFF, zip, audit JSON
src/gemsdoe40/           ← metric, Euler, KDE, uniqueness, holdout
scripts/run_euler_depth_cluster.py
tests/                   ← metric worked example, synthetic Euler, format contract
```

## Licence / data

Competition rasters remain under the DrivenData / USGS / INGENIOUS
licences. GeoDAWN magnetics are CC0 (USGS DOI 10.5066/P93LGLVQ). Code in
this repository is for the GEMS Prize and is not a vendor product.
