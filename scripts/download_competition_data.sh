#!/usr/bin/env bash
# Assemble the official competition rasters.
# Preferred: concatenate the GitHub-hosted 5GEMSDOE data-bridge parts
# (sha256-pinned to the official training_features.tif).
# Fallback: DrivenData login download into data/.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="$ROOT/data"
mkdir -p "$DATA/bridge"

echo "This script does not log into DrivenData."
echo "Official tab (login-walled): https://www.drivendata.org/competitions/306/competition-doe-gems/data/"
echo
echo "If data/bridge/gems-geodawn-numerical-features.tif.part-* exist, assembling..."

PARTS=(
  "$DATA/bridge/gems-geodawn-numerical-features.tif.part-000"
  "$DATA/bridge/gems-geodawn-numerical-features.tif.part-001"
  "$DATA/bridge/gems-geodawn-numerical-features.tif.part-002"
  "$DATA/bridge/gems-geodawn-numerical-features.tif.part-003"
  "$DATA/bridge/gems-geodawn-numerical-features.tif.part-004"
)
if [[ -f "${PARTS[0]}" ]]; then
  cat "${PARTS[@]}" > "$DATA/training_features.tif"
  echo "assembled $DATA/training_features.tif ($(wc -c < "$DATA/training_features.tif") bytes)"
  python3 - << 'PY'
import hashlib, pathlib
p = pathlib.Path("data/training_features.tif")
h = hashlib.sha256(p.read_bytes()).hexdigest()
exp = "4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5"
print("sha256", h)
print("match", h == exp)
if h != exp:
    raise SystemExit("hash mismatch — refuse to proceed")
PY
else
  echo "parts not present. Download training_features.tif from DrivenData into data/"
  echo "or fetch the 5GEMSDOE data-bridge parts (see data/bridge/manifest.json)."
fi
