# Session 3 — H11 frozen pre-screen: NEGATIVE-AT-PRESCREEN (2026-10-06)

Status: research record only. No submission slot was used. No organizer score is claimed.

## What was tested

**H11** (registered rank 4 in the session-2 slate): seismicity-corridor × shallow-Euler
intersection — cells simultaneously near recorded seismicity (`deq_n100a15` band 10,
`ieq_n100a15` band 16) and under shallow depth-consistent Euler support (frozen H4 cloud).
This was the last locally testable registered candidate (H12 remains egress-blocked).

Same fail-fast pattern: protocol registered before scoring (preregistration addendum 4,
pinned at sha256 `d6117108c59f18371f25a428b608193d98abc0b03fd13f1c018d16b9c76a2a99`;
`scripts/prescore_h11_bands.py` fails closed on document drift, and additionally verifies
the frozen cloud sha `6bed30b2…`).

## Protocol (frozen before results)

- σ = 2 px pre-smoothing for derivatives; six frozen indicators: percentile proximity to
  seismicity (−band 10, −band 16), gradient magnitudes of both distance bands, the rank
  product of the two proximities, and **`prox_deq_x_euler`** — the registered intersection
  signature itself (proximity rank × Euler-support KDE rank, KDE built with the exact
  H8-family construction: bilinear splat + σ = 2, truncate = 4 Gaussian).
- Percentile emission on the template-defined 5,167,373-cell footprint; frozen SGMC proxy
  (26d142…); the audits' exact `mass_matched_controls`.
- Decision rule: advance only if the best indicator beats BOTH 0.039354 and 2 × its own
  best random control.

## Result

| indicator | proxy score | random control | gradient-topK |
|---|---|---|---|
| prox_deq_x_euler | **0.0748** | 0.1156 | 0.0800 |
| HG_deq | 0.0702 | 0.1158 | 0.0800 |
| HG_ieq | 0.0689 | 0.1158 | 0.0800 |
| prox_both | 0.0639 | 0.1868 | 0.0832 |
| prox_deq | 0.0637 | 0.1868 | 0.0832 |
| prox_ieq | 0.0631 | 0.1868 | 0.0832 |

**No indicator reaches its own random control.** The registered intersection signature is
the strongest of the six, yet still far below random. Decision:
**NEGATIVE-AT-PRESCREEN**. H11 is closed. Full receipt: `docs/data/h11-prescore.json`.

## Interpretation (registered caveats apply)

1. Recorded-seismicity proximity (even intersected with frozen Euler support) does not
   place mass on the off-catalogue SGMC traces better than random dots at equal budget.
   The proximity fields' own random controls run high (0.1868) because their emissions
   concentrate where the proxy also has mass, but the fields never beat them.
2. With H9, H10 and H11 all closed at pre-screen, **every locally testable hypothesis in
   the session-2 register has now been evaluated**: two implemented and held at the frozen
   gates (H8, H13), three closed negative at pre-screen (H9, H10, H11), one egress-blocked
   (H12). The register is exhausted without finding placement skill.
3. The stable empirical fact across all five pre-screens/audits: the frozen proxy rewards
   diffuse mass (random ≈ 0.116 at this budget), and no local signal — Euler, subsurface
   bands, geodetic strain, or seismicity proximity — beats it. What separates live
   leaderboard scores is alignment with hidden traces that cannot be reconstructed locally.
4. All inputs are public-mirror rasters and a historical SGMC-derived proxy; results are
   statements about these instruments, not the hidden test set.

## Consequences

- Hypothesis register: H11 closed NEGATIVE-AT-PRESCREEN; register exhausted except the
  blocked H12.
- H8, H13, H8-ASA (concurrent session) all remain HOLD — DO NOT SUBMIT. No weekly slot
  used; no organizer score claimed.
- The binding prerequisite for any future slot spend is unchanged and now better
  evidenced: a registered, source-independent promotion instrument must exist before any
  further emission work. New ideas for the next session should be registered against that
  gap (e.g., external authenticated data, or a non-proxy validation design), not as
  further band re-combinations on the same frozen proxy.
