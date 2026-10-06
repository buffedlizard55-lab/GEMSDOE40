# Prior-raster inventory refresh — 2026-10-06

## Scope and result

The refreshed public-owner scan covered **55 GEMSDOE repositories at their recorded default-branch heads**. All 55 recursive Git trees were read without truncation. It found **501 TIFF paths**; every current-head path received an explicit include/exclude classification, with **zero unclassified paths**. The scan is an audit of public repository artifacts, not a complete census of private or uncommitted submissions.

| Current-head TIFF path classification | Paths | Decision |
|---|---:|---|
| Prior prediction/submission/research outputs | 358 | Include |
| Explicit `submission.tif` / `submission_conformant.tif` paths previously marked “likely” | 48 | Manually reviewed from repository/path and filename context; include as prior prediction/submission outputs |
| Prior prediction outputs in docs root | 8 | Include |
| Model probability outputs (variants) | 4 | Include; retain variant identity |
| Proxy-circular research/prediction outputs | 4 | Include and flag for interpretation |
| Input/template/label/proxy/source rasters | 76 | Exclude from prediction corpus |
| Input/template/label/proxy/source fixtures | 3 | Exclude from prediction corpus |
| **Total** | **501** | **0 pending / unclassified** |

The 48 formerly “likely” paths were all named exactly `submission.tif` or `submission_conformant.tif` under archived baseline/run/union evidence directories; they are now directly classified as prior submission/prediction outputs. This resolves the pending row-level path review without treating every TIFF in the repository as a prediction.

The refreshed corpus contains **306 unique Git blobs** across **541 pinned repository/ref/path locations**, totaling **487,777,241 bytes** over the unique blobs. The previous inventory had 279 blobs; this refresh adds **27 previously unseen unique prediction-output blobs**. The 82 current-head prediction locations added or confirmed by the refresh include repeated paths and are not 82 new unique artifacts. The full pinned identities, paths, classifications, SHA-256 values, and repository heads are in [`prior_raster_inventory-20261006.json`](../data/prior_raster_inventory-20261006.json).

## Byte/cache verification

All **306/306 unique TIFF blobs** were fetched into the external temporary cache and verified by `scripts/fetch_prior_cache.py` against their recorded Git blob SHA-1 and SHA-256. Cache: `/tmp/gemsdoe40-prior-cache` (not committed). The 27 new blobs were additionally byte-checked during the inventory refresh itself. GEMSDOE40 is one of the 55 scanned repositories, but the refreshed inventory pins its recorded public default-branch head from before H7 was built. The unpushed H7 working-branch candidate is therefore not in the prior corpus.

## Limitations and interpretation

- Repository paths are not proof that a raster was uploaded to DrivenData, scored, selected, or associated with a public leaderboard account. The scan must not be described as an official submission ledger.
- The corpus intentionally excludes feature, label, sample/template, raw proxy, source, and test fixture TIFFs. Some research outputs are not portal-format-valid or were never scored; they remain useful for an artifact-level uniqueness audit.
- Four identified proxy-circular output paths are retained rather than dropped. They must be visibly flagged in any SGMC-proxy score comparison; their apparent performance is not independent evidence of hidden-set performance.
- A raster blob is fetched once even if the same bytes occur at many pinned locations. The inventory retains the separate locations, while pixel comparisons operate on unique blobs.
- Completeness is bounded to public owner repositories/default-branch trees at the recorded commits and to paths classified by the checked-in rules. It cannot establish uniqueness against private repositories or unseen future outputs.

## Next use

Use this exact 306-blob corpus for H7 canonical-pixel/raw-raster uniqueness comparisons and for the frozen same-grid holdout incumbent scan. Require a complete cache before declaring uniqueness. Keep all current/historical proxy versions pinned separately; neither owner mirror is authenticated organizer truth.
