#!/usr/bin/env python3
"""Compatibility entry point: rebuild the current site, then assert the retained H4 identity.

The H4 pages were retired on 2026-10-06, first in favour of the session-2/3 H13 / H8
trace-locked / H8-ASA pages, and now in favour of the H8 lineament release. CI, the Pages
runner and older documentation still call this path, so it delegates to the single current
generator (``scripts/build_h8_site.py``) and fails closed if the retained H4 artifact bytes
were ever modified. The retired session-2/3 generator is kept for audit in
``scripts/retired/build_h13_site_session2.py``.
"""
from __future__ import annotations
import hashlib
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import build_h8_site  # noqa: E402

H4 = ROOT / "docs/downloads/gemsdoe40-h4-contact-offset-depthkde-20261006-fab9f6619c02.tif"
H4_SHA256 = "ee73ffd76fabbaa1a2e77e17f57947a7db858916d713801e0c3502e49f6acabb"


def assert_retained_h4() -> None:
    """Fail closed if the retained H4 release bytes were ever modified."""
    actual = hashlib.sha256(H4.read_bytes()).hexdigest()
    if actual != H4_SHA256:
        raise SystemExit(f"retained H4 artifact changed: {actual} != {H4_SHA256}")


if __name__ == "__main__":
    assert_retained_h4()
    raise SystemExit(build_h8_site.main())
