#!/usr/bin/env bash
# Fetch the 17 prior/incumbent rasters that the novelty audit and the augmentation are measured
# against, from the repositories that published them, into data/prior/ (git-ignored).
#
# These are *other repositories' outputs*, not organizer data.  They are fetched only to (a) prove
# this repository's artifacts are not near-duplicates of anything already submitted and (b) build
# the measured augmentation baseline.  No prior raster is uploaded anywhere by this repository.
#
# Requires: gh (authenticated).  Run: bash scripts/fetch_priors.sh
set -u
cd "$(dirname "$0")/.."
mkdir -p data/prior

fetch() {  # repo path dest
  gh api "/repos/buffedlizard55-lab/$1/contents/$2" -H "Accept: application/vnd.github.raw" > "data/prior/$3"
  echo "$(sha256sum "data/prior/$3" | cut -c1-16) $(stat -c%s "data/prior/$3") $3"
}

fetch GEMSDOE24 data/external/derived_sgmc_faults_100m_u8.tif sgmc_faults_u8.tif
fetch GEMSDOE32 docs/downloads/gemsdoe32-h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros.tif top_02778_h33b2_zeros.tif
fetch GEMSDOE32 docs/downloads/gemsdoe32-h33-h33-2-b2-20261004T220000Z-e5eb6e7e-nan.tif   top_02778_h33b2_nan.tif
fetch GEMSDOE28 docs/downloads/gems28-h27-4-r1-solo-d2-8-20261003-8acb75e1f2cc-nan.tif        base_02708_h274r1_d28_nan.tif
fetch GEMSDOE24 docs/downloads/gems24-h25-1-dotted-h19-5-d1-5-20261002-989f59505db1-nan.tif   prior_d15_h195_nan.tif
fetch GEMSDOE24 docs/downloads/gems24-h25-1-dotted-h19-5-d2-8-20261002-e56ea318af89-nan.tif   prior_d28_h195_nan.tif
fetch 19GEMSDOE  docs/downloads/gems19-h19-5-powerlaw-budget-multiline-corroborated-20260930-e27054cf-nan.tif prior_h195_solid_nan.tif
fetch 16GEMSDOE  docs/downloads/gems16-h16-1-topo-geophys-baseline-ridges-20260930-df20f65e-nan.tif   prior_h16_1_nan.tif
fetch GEMSDOE10  docs/downloads/gems10-h25-ctx-ridge-20260927T232947704150Z-6452ae1d00.tif      prior_h25ctx_nan.tif
fetch GEMSDOE10  docs/downloads/gems10-h28-dotted-ridge-20260928T020256236880Z-6452ae1d00.tif   prior_h28dotted_nan.tif
fetch 12GEMSDOE  docs/downloads/12GEMSDOE_r7-nms3-dem10-scarp_0c9199f14e62_allfinite.tif       prior_r7scarp_allfinite.tif
fetch 13GEMSDOE  docs/downloads/13gems_20261001_r13-lattice-s5_v2_nan-outside.tif             prior_lattice_nan.tif
fetch GEMSDOE24  inputs/calibration/8GEMSDOE_Hedge-v2_submission.tif                           prior_hedgev2.tif
fetch GEMSDOE28  docs/downloads/gems27-topo-gap-closure-t-v2-on-d1-5-20261002-5512495c6bd1-nan.tif prior_tv2_d15_nan.tif
fetch GEMSDOE28  docs/downloads/gems28-h32-1-prethin-tip-euler-d2-8-20261003-31e35eee884e-nan.tif  prior_h32_1_tipeuler_nan.tif
fetch GEMSDOE24  inputs/calibration/gemsdoe9-PLACEHOLDER-2314b599.tif                          prior_placeh_nan.tif
fetch GEMSDOE24  inputs/calibration/gemsdoe-ens12-adopted-7f00890a.tif                          prior_ens12_nan.tif

echo "--- data/prior/ has $(ls -1 data/prior/*.tif | wc -l) rasters ---"
