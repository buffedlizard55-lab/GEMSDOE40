# Session 3 — H9 frozen pre-screen: NEGATIVE-AT-PRESCREEN (2026-10-06)

Status: research record only. No submission slot was used. No organizer score is claimed.

## What was tested

**H9** (registered rank 2 in the session-2 slate): blind basement-flexure faults, using the
never-consumed feature bands `depth_to_base_surf` (15), `cond_surf` (17),
`iso_grav_anom_hg` (18), with `det_elev` (12) as registered companion.

Instead of building a full emission pipeline first, session 3 ran a **frozen pre-screen**:
one registered protocol, one run, one decision. The protocol was appended to
`docs/research/h8-preregistration-20261006.md` as *addendum 2* and pinned at sha256
`95ef16f16cafc0926ae02339676a76f21a826d19a2b0e8f23f279e6bdb47c2ee` **before** any H9
scoring; `scripts/prescore_h9_bands.py` refuses to run if that document changes.

## Protocol (frozen before results)

- Inputs: pinned `data/training_features.tif` (sha256 `4371c82e…`), pinned
  `data/sample_submission.tif` template (sha256 `2176d08e…`) defining the 5,167,373-cell
  competition footprint, frozen SGMC proxy `derived_sgmc_faults_100m_u8.tif` (26d142…,
  hash-enforced inside `read_proxy_truth`).
- Smoothing: single Gaussian σ = 2 px (200 m) per raw band before any derivative.
- Indicators (exactly six): HG of bands 15, 17, 18, 12; percentile-rank product of
  HG(15) × HG(17) ("flexure_product"); |Laplacian| of band 15 ("step_depth_base").
- Emission: each indicator percentile-normalized to [0, 1] over finite in-footprint cells
  (measurement instrument only — never a submission candidate).
- Controls: the audit's exact `mass_matched_controls` (seeds 40/41/42 random dots at the
  value-mass budget + TMI gradient-top-K), scored with `score_array_on_proxy`.
- Decision rule: advance H9 only if the best indicator beats **both** 0.039354
  (2 × H13's session-2 G3 score 0.019677) **and** 2 × its own best random control.

## Result

| indicator | proxy score | random control | gradient-topK |
|---|---|---|---|
| HG_det_elev | **0.1048** | 0.1158 | 0.0800 |
| HG_cond | 0.0710 | 0.1158 | 0.0800 |
| HG_grav_hg | 0.0680 | 0.1158 | 0.0800 |
| step_depth_base | 0.0650 | 0.1156 | 0.0800 |
| flexure_product | 0.0619 | 0.1158 | 0.0800 |
| HG_depth_base | 0.0579 | 0.1158 | 0.0800 |

Every indicator passes the absolute bar (all > 0.039354), but **no indicator reaches its
own random control, let alone 2 × it**. Decision: **NEGATIVE-AT-PRESCREEN**. H9 is closed;
no implementation compute will be spent on it. Full receipt with per-indicator controls,
masses, hashes and timestamps: `docs/data/h9-prescore.json`.

## Interpretation (registered caveats apply)

1. The three subsurface bands central to H9 (15/17/18) are the **weakest** indicators;
   only the detrended-elevation gradient (band 12) comes near random, still below it.
2. At dense emission mass (value budget 2,577,859 px) random dots score 0.1158 on the
   proxy — consistent with the session-2 finding that the SGMC proxy rewards diffuse mass
   distribution, while what separates live scores is placement on unknowable hidden traces.
3. The proxy's incumbent ceiling (0.8359, SGMC-matching fields) remains circular and was
   not used as any bar here.
4. Bands 15/17/18 come from public mirrors, not organizer-authenticated downloads; the
   negative result is a statement about these mirrors and this frozen proxy.

## Process note (disclosed)

The first pre-screen attempt derived the footprint from feature-band finiteness and so
emitted over the full 12,279,160-px rectangle instead of the sample template's
5,167,373-cell competition footprint. That run's controls were not mass-matched to the
competition footprint and are **void**. The bug was fixed (footprint now read from the
pinned `sample_submission.tif`) and the corrected run above is the only recorded result.

## Consequences

- Hypothesis register: H9 closed NEGATIVE-AT-PRESCREEN. Next registered candidates:
  **H10** (geodetic strain-rate lineaments, bands 4/7/8 — cheap, same pre-screen pattern
  applies), H11 (seismic corridor × Euler), H12 (ComCat; egress-blocked here).
- No weekly slot was used; H8 and H13 remain HOLD — DO NOT SUBMIT; the repo's
  source-independent promotion protocol is still unbuilt and remains the binding
  prerequisite for any future slot spend.
- Append-only discipline preserved: addendum 2 was appended to the preregistration before
  scoring; the H8 and H13 validation receipts were re-audited and now disclose the
  amendment (generation-time hashes reproduced: `6efaaf64…` for H8, `5fc48198…` for H13).
