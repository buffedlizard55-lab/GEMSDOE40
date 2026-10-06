"""Prevent accidentally publishing archived model variants as current advice."""
from __future__ import annotations
import os
import sys


def require_legacy_opt_in() -> None:
    message = (
        "ARCHIVED RESEARCH ENTRY POINT: not the current candidate or promotion protocol. "
        "Use scripts/run_contact_euler.py followed by scripts/audit_contact_candidate.py. "
        "For deliberate historical reproduction only, set GEMSDOE_ALLOW_LEGACY_RESEARCH=1. "
        "Legacy scores/uniqueness checks are not approval to submit."
    )
    if os.environ.get("GEMSDOE_ALLOW_LEGACY_RESEARCH") != "1":
        raise SystemExit(message)
    print(message, file=sys.stderr)
