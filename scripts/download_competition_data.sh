#!/usr/bin/env bash
# Autonomous restoration from public, immutable, hash-verified mirrors.
# This does not authenticate to DrivenData or claim organizer provenance.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
if [[ -n "${PYTHON:-}" ]]; then
  PY="$PYTHON"
elif [[ -x "$ROOT/.venv/bin/python" ]]; then
  PY="$ROOT/.venv/bin/python"
else
  PY=python3
fi
"$PY" scripts/acquire_data.py "$@"
