#!/usr/bin/env bash
# Compatibility entry point. The old mutable/append-in-place downloader was unsafe.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
exec bash "$ROOT/scripts/download_competition_data.sh" "$@"
