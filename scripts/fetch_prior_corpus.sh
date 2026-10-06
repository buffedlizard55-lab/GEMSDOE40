#!/usr/bin/env bash
# Fetch the public prior-submission corpus used for calibration and novelty audit.
#
# Each entry is a GeoTIFF published by a sibling GEMSDOE repository.  The files
# are NOT committed (data/prior/ is gitignored); this script restores them from
# pinned repository paths using a blobless sparse clone, so no file is ever
# copied into a submission.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEST="${1:-$ROOT/data/prior}"

# repo:path-with-sparse-dirs
REPOS=(
  "GEMSDOE32:docs/downloads"
  "GEMSDOE28:docs/downloads"
  "GEMSDOE33:docs/downloads"
  "GEMSDOE24:docs/downloads data/external"
  "GEMSDOE27:docs/downloads"
  "GEMSDOE25:docs/downloads"
  "GEMSDOE29:docs/downloads"
  "GEMSDOE30:docs/downloads"
  "GEMSDOE31:docs/downloads"
  "19GEMSDOE:docs/downloads"
  "16GEMSDOE:docs/downloads"
  "GEMSDOE10:docs/downloads"
  "13GEMSDOE:data/external_sgmc"
)

mkdir -p "$DEST"
for entry in "${REPOS[@]}"; do
  repo="${entry%%:*}"; dirs="${entry#*:}"
  target="$DEST/$repo"
  if [ ! -d "$target/.git" ]; then
    echo "cloning $repo (blobless, shallow)"
    git clone --filter=blob:none --no-checkout --depth=1 \
      "https://github.com/buffedlizard55-lab/$repo" "$target" >/dev/null 2>&1 || {
        echo "  FAILED to clone $repo"; continue; }
    git -C "$target" sparse-checkout init --cone >/dev/null 2>&1
  fi
  # shellcheck disable=SC2086
  git -C "$target" sparse-checkout set $dirs >/dev/null 2>&1
  git -C "$target" checkout >/dev/null 2>&1 && echo "  $repo: $(du -sh "$target" | cut -f1)"
done
echo "prior corpus under $DEST"
