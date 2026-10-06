#!/usr/bin/env python3
"""Hard-link the pinned prior corpus into the layout ``gemsdoe40.corpus`` expects.

``scripts/fetch_prior_corpus.sh`` clones each sibling repository (blobless, sparse)
under ``data/prior/<repo>/``.  The calibration and novelty-audit code reads
``data/prior/<key>/<filename>`` instead, so this step maps the two.  Hard links
mean no extra bytes are written; a re-run is idempotent and verifies every source
file exists before linking.

Usage: python3 scripts/link_prior_corpus.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe40.corpus import PRIORS, PRIOR_DIR, local_path  # noqa: E402


def main() -> int:
    linked = present = missing = 0
    for p in PRIORS:
        dst = local_path(p)
        if dst.exists():
            present += 1
            continue
        src = PRIOR_DIR / p.repo / p.path
        if not src.exists():
            print(f"MISSING SOURCE  {p.key:34s} {p.repo}/{p.path}")
            missing += 1
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        os.link(src, dst)
        linked += 1
    print(f"priors: {len(PRIORS)} · already present {present} · linked {linked} · missing {missing}")
    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
