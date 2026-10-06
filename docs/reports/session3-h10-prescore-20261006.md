# Session 3 — H10 frozen pre-screen: NEGATIVE-AT-PRESCREEN (2026-10-06)

Status: research record only. No submission slot was used. No organizer score is claimed.

## What was tested

**H10** (registered rank 3 in the session-2 slate): geodetic strain-rate lineaments from
`geod_shearrate` (7), `geod_dilaterate` (8), `geod_2ndinv` (4) — the registered signature
being oriented structure-tensor lineament extraction, with relay segments between catalogue
traces.

Same fail-fast pattern as the H9 pre-screen: one frozen protocol, one run, one decision,
registered before scoring (preregistration addendum 3, pinned at sha256
`d495fc0f91509ea41120514ab0fd2b10aac4215f2d0c55ac5abb36533569bea3`;
`scripts/prescore_h10_bands.py` fails closed on document drift).

## Protocol (frozen before results)

- Inputs: pinned feature raster `4371c82e…`, pinned template `2176d08e…` defining the
  5,167,373-cell competition footprint (the H9 footprint bug was avoided from the start),
  frozen SGMC proxy 26d142… hash-enforced.
- σ = 2 px Gaussian before derivatives; six frozen indicators: HG of bands 7, 8, 4;
  band 4 used directly (`inv_direct`); structure-tensor oriented linearity of band 7
  (coherence × √λ1, tensor smoothed σ = 2); percentile-rank product HG(7) × HG(8).
- Emission percentile-normalized to [0,1]; scored with `score_array_on_proxy` plus the
  audits' exact `mass_matched_controls`.
- Decision rule: advance only if the best indicator beats BOTH 0.039354 and 2 × its own
  best random control.

## Result

| indicator | proxy score | random control | gradient-topK |
|---|---|---|---|
| HG_dilate | **0.0701** | 0.1158 | 0.0800 |
| HG_2ndinv | 0.0690 | 0.1158 | 0.0800 |
| inv_direct | 0.0675 | 0.1156 | 0.0800 |
| shear_x_dilate | 0.0668 | 0.1158 | 0.0800 |
| lineament_shear | 0.0653 | 0.1156 | 0.0800 |
| HG_shear | 0.0647 | 0.1158 | 0.0800 |

**No indicator reaches its own random control**, let alone 2 × it. Decision:
**NEGATIVE-AT-PRESCREEN**. H10 is closed; no implementation compute will be spent.
Full receipt: `docs/data/h10-prescore.json`.

## Interpretation (registered caveats apply)

1. The registered oriented-lineament detector (`lineament_shear`) is second-weakest of the
   six; plain gradient magnitudes of the strain bands add nothing either. Geodetic
   strain-rate on these mirrors does not localize off-catalogue SGMC fault traces beyond
   what random mass placement achieves.
2. The dense-mass random level (≈0.1156–0.1158) reproduced exactly across H9 and H10 —
   the proxy's reward for diffuse mass is stable and candidate-independent.
3. With H9 and H10 both closed at pre-screen, every registered "new band family" lever
   (subsurface + geodetic) has now failed the proxy placement test. Remaining registered
   candidates: H11 (seismic corridor × Euler intersection — uses spatial proximity to
   recorded seismicity rather than band gradients) and H12 (ComCat; egress-blocked here).
4. All inputs are public-mirror rasters, not organizer-authenticated; the negative result
   is a statement about these mirrors and this frozen proxy.

## Consequences

- Hypothesis register: H10 closed NEGATIVE-AT-PRESCREEN.
- H8 (trace-locked) and H13 (crest-binary) remain HOLD — DO NOT SUBMIT; H8-ASA (concurrent
  session) likewise HOLD. No weekly slot was used; no organizer score is claimed.
- The binding prerequisite for any future slot spend remains a registered,
  source-independent promotion protocol; none exists yet.
