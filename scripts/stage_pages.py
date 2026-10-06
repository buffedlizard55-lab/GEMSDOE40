#!/usr/bin/env python3
"""Stage only public docs; preserve both /GEMSDOE40/docs/ and old root URLs."""
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]


def stage() -> Path:
    target = ROOT / "_site"
    # Fixed generated directory only: never a user-selected or repository root.
    if target.exists(): shutil.rmtree(target)
    target.mkdir()
    shutil.copytree(ROOT / "docs", target / "docs")
    shutil.copytree(ROOT / "docs", target, dirs_exist_ok=True)
    shutil.copyfile(ROOT / "index.html", target / "index.html")
    (target / ".nojekyll").touch()
    for required in (target / "index.html", target / "docs/index.html", target / "docs/data/current-candidate.json"):
        if not required.is_file(): raise RuntimeError(f"missing Pages route: {required}")
    return target


if __name__ == "__main__":
    print(f"Staged public content at {stage()}; no .git, credentials, input rasters or prior cache copied")
