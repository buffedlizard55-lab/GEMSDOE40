#!/usr/bin/env python3
"""Deprecated name kept for CI: rebuild the current site, then assert the retained H4 identity.

The H4 pages were retired on 2026-10-06 when the H8 lineament candidate replaced them as the
current artifact.  H4 is still published as *retained history* and must never be altered, so
this wrapper calls ``scripts/build_h8_site.py`` (the single deterministic builder) and then
re-checks the H4 release bytes by SHA-256.
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
