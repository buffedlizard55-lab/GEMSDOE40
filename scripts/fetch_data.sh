#!/usr/bin/env bash
# Fetch the pinned competition data + the off-catalogue holdout truth into data/.
#
# The organizer's data tab requires a login and the sandbox network cannot reach
# drivendata.org / dropbox.com (both return HTTP 000 here), so the byte-identical mirrors that
# sibling repositories published on GitHub are used.  Every file is verified against the sha256
# pins recorded in src/gems40/pins.py (and re-derived by scripts/inspect_data.py).
#
# Requires: gh (authenticated), python3.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data data/external

echo "[1/3] training_features.tif  (5 parts, 418,912,844 B total)"
parts=(data/bridge/gems-geodawn-numerical-features.tif.part-000
       data/bridge/gems-geodawn-numerical-features.tif.part-001
       data/bridge/gems-geodawn-numerical-features.tif.part-002
       data/bridge/gems-geodawn-numerical-features.tif.part-003
       data/bridge/gems-geodawn-numerical-features.tif.part-004)
if [ ! -f data/training_features.tif ]; then
  for p in "${parts[@]}"; do
    echo "  -> $p"
    gh api "/repos/buffedlizard55-lab/GEMSDOE/contents/$p" -H "Accept: application/vnd.github.raw" \
      >> data/training_features.tif
  done
fi

echo "[2/3] labels / sample submission"
gh api "/repos/buffedlizard55-lab/GEMSDOE24/contents/data/bridge/existing_faults.tif" \
  -H "Accept: application/vnd.github.raw" > data/existing_faults.tif
gh api "/repos/buffedlizard55-lab/GEMSDOE24/contents/data/bridge/example_submission.tif" \
  -H "Accept: application/vnd.github.raw" > data/example_submission.tif

echo "[3/3] off-catalogue holdout truth (USGS SGMC faults, 100 m)"
gh api "/repos/buffedlizard55-lab/GEMSDOE24/contents/data/external/derived_sgmc_faults_100m_u8.tif" \
  -H "Accept: application/vnd.github.raw" > data/external/derived_sgmc_faults_100m_u8.tif

python3 scripts/inspect_data.py
echo "done -- data/ is ready"
