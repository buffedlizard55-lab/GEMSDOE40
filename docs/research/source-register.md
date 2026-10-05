# Evidence and source register — GEMSDOE40

**Last checked:** 2026-10-05 UTC. A link to an official source does not authenticate a local mirror or prove that a layer is fit for a particular analysis.

## Official competition sources

| Source | What it establishes | Limit / handling |
|---|---|---|
| [DrivenData task and metric](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/) | Task is geologic fault prediction in the GeoDAWN region, with a distance-weighted Tversky metric. | Does not provide hidden labels or identify which file earned an account-level leaderboard score. |
| [DOE/NLR official rules PDF](https://docs.nlr.gov/docs/fy26osti/96647.pdf) | Official competition/submission rules, including output requirements and schedule. | Portal authentication/enrollment is required for actual submission. Rules have timing language that should be checked before any future upload. |
| [DrivenData staff clarification](https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516/2) | The supplied known-fault mask is removed pixel-exactly; nearby off-mask pixels remain scoreable. | This is a staff forum response, not a substitute for the full rules or a guarantee about final-round labels. |
| [Official public leaderboard](https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/) | Snapshot read on 2026-10-05: top four 0.3262, 0.3222, 0.3220, 0.3195; `extradr19` has a 0.2778 account-level row. | The row does not expose a TIFF-level receipt. Do not attribute 0.2778 to H33-2-B2; the GEMSDOE32 record says 0.2747 projected and UNSCORED. Dated local copy: [`feed.json`](../data/feed.json). |
| [DrivenData data tab](https://www.drivendata.org/competitions/306/competition-doe-gems/data/) | Intended authenticated route for official feature/label/template rasters and other challenge data. | Unauthenticated sandbox access redirected to login. Local project inputs arrived through a public sibling-repository bridge and are hash-pinned in the experiment reports; they are not represented as an authenticated download. |

## Scientific / public data sources

| Source | Relevance | Limit / handling |
|---|---|---|
| [USGS SGMC metadata](https://mrdata.usgs.gov/geology/state/USGS_SGMC_Metadata.xml), [USGS Data Series 1052](https://pubs.usgs.gov/ds/1052/) | Official metadata and provenance for the vector State Geologic Map Compilation; includes mapped structural features. | The evaluation raster is an owner-derived mirror. It was not independently rebuilt from the USGS vectors. Its exact provenance and rasterization are unverified; it is only a proxy instrument, never hidden challenge truth. |
| [USGS ComCat FDSN Event Web Service](https://earthquake.usgs.gov/fdsnws/event/1/) | Official source with event details and associated moment-tensor/focal-mechanism products. | The locked sample-extent query returned 34,716 event records advertising a `focal-mechanism` product (238 at magnitude ≥4); one event detail exposed finite nodal planes. The usable-plane count and distribution across 24 blocks remain unverified, so H3 is still blocked. Direct sandbox HTTPS failed; no ComCat data were used to generate H1/H2/H2-B. See [`h3_source_audit.json`](../data/h3_source_audit.json). |
| Reid & Thurston (2014), [“The structural index in gravity and magnetic interpretation: Errors, uses, and abuses”](https://www.reid-geophys.co.uk/wp-content/uploads/2017/11/Reid-Thurston-2014.pdf) | Structural-index guidance. SI=0 is cited for an infinite magnetic contact; gravity SI=−1 is cited for an idealized finite contact/fault. | H1's SI=0 on the first vertical derivative of gravity is explicitly an approximation, not exact finite-contact physics. |
| [GEMSDOE28 Euler research](https://buffedlizard55-lab.github.io/GEMSDOE28/) | Sibling-project prior art: TMI-only SI=0 Euler clouds and later heat-flow evidence. | Not official and not independent validation. GEMSDOE40's novelty claim is narrow: cross-field/depth consensus or multi-height persistence, not Euler itself. |

## Local artifact identity and irregularities

- Training features SHA-256: `4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5`.
- Example template SHA-256: `2176d08e485aa2cd2860ce8df539db4faf4d76163b38a4dd8c30a40454d35cbc`.
- Existing-label mask SHA-256: `7ba308ccdc4418b31a178f4f1ef21aaa6e152e4028f2f6f64b01f7eb25ae4093`.
- Owner-derived SGMC proxy SHA-256: `26d142c4c93282cd94f6950ab96f22aeff59fbbea523d43d662e76fa1b161b5c`.
- The mirrored example template contains positive known-label pixels even though the official page describes an all-zero example. This project uses the example's grid/footprint and the separate labels raster; it does not treat the sample's pixel values as target truth. This discrepancy is preserved as an irregularity, not silently normalized.
- H1's magnetic vertical-gradient check passed (`r=0.9959886`, scale ratio `0.988864`); the supplied gravity vertical-gradient feature did not agree with the raw-gravity Fourier derivative (`r=0.17839`, std 1.628 versus 0.00266 computed), so H1 derived gravity gradients from raw `iso_grav_anom`. These are operator diagnostics, not proof of correct geologic source physics.
- H1 stopped with zero gravity solutions; strict H2 stopped below its fixed support budget. H2-B's separate output rule and negative proxy holdout are fully documented. No post-result threshold change is hidden.

## Prior-raster comparator caveat

The all-prior inventory contains 270 unique TIFF Git blobs; 268 are on the exact sample grid, and 166 are format-eligible. The best format-eligible prior proxy score, 0.8359066541, belongs to the pinned `GEMSDOE3` artifact named `gapfinder-v2-sgmc-gap`. Its originating project describes the output as an SGMC-gap-only raster. Since the local proxy is derived from SGMC, this prior is circular on the proxy. The frozen gate retains it exactly as registered; this makes the candidate's gate failure unambiguous but does **not** make 0.8359 an independent estimate of hidden-label generalization.

The circularity note was appended after the initial score computation; no score values or gate logic were changed. The same provenance warning is stored in [`validation.json`](../data/validation.json).
