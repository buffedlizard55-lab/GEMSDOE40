#!/usr/bin/env python3
"""Build the complete owner-reported score register and pinned project-page review.

Reads public GitHub page SOURCE at immutable commits. This is not represented
as an HTTP uptime audit or an organizer verification of user-reported scores.
"""
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]


class Text(HTMLParser):
    def __init__(self):
        super().__init__(); self.parts = []; self.links = []; self.skip = 0
    def handle_starttag(self, tag, attrs):
        if tag in ("style", "script"): self.skip += 1
        if tag == "a":
            a = dict(attrs)
            if a.get("href"): self.links.append(a["href"])
    def handle_endtag(self, tag):
        if tag in ("style", "script"): self.skip = max(0, self.skip-1)
    def handle_data(self, data):
        if not self.skip and data.strip(): self.parts.append(data.strip())


def parse_prompt(text: str) -> list[dict]:
    prefix = text.split("WE NEED TO STUDY, ANALYZE, AND UNDERSTAND THE HIGHEST SCORE")[0]
    groups, current = [], None
    for line in prefix.splitlines():
        match = re.match(r"\[(https://buffedlizard55-lab\.github\.io/([^/]+)/[^\]]*)\]", line)
        if match:
            current = {"repo": match[2], "site": match[1], "submissions": []}
            groups.append(current); continue
        alias = re.fullmatch(r"(4[0-4])GEMSDOE", line.strip())
        if alias:
            repo = "GEMSDOE" + alias[1]
            current = {"repo": repo, "site": f"https://buffedlizard55-lab.github.io/{repo}/",
                       "site_note": "URL inferred from the bare project alias in the brief", "submissions": []}
            groups.append(current); continue
        result = re.fullmatch(r"([A-Za-z0-9_][A-Za-z0-9_.-]*):\s*(0\.\d+)?\s*", line)
        if result and current:
            current["submissions"].append({"name": result[1], "score": float(result[2]) if result[2] else None,
                                           "evidence": "owner-reported" if result[2] else "no score supplied"})
    return groups


def main():
    projects = parse_prompt((ROOT / "docs/user-prompt-20261006.md").read_text())
    inventory = json.loads((ROOT / "docs/data/prior-inventory-20261006.json").read_text())
    snapshots = {s["repo"]: s["commit"] for s in inventory["repository_snapshots"]}
    def review(project):
        repo, ref = project["repo"], snapshots[project["repo"]]
        attempts = []
        for path in ("docs/index.html", "index.html", "README.md"):
            command = ["gh", "api", "-H", "Accept: application/vnd.github.raw", f"repos/buffedlizard55-lab/{repo}/contents/{path}?ref={ref}"]
            response = subprocess.run(command, capture_output=True, timeout=60)
            if response.returncode:
                attempts.append({"path": path, "status": "not returned by GitHub"}); continue
            html = response.stdout.decode("utf-8", errors="strict")
            parser = Text(); parser.feed(html)
            text = " ".join(parser.parts) if path.endswith(".html") else html
            title_match = re.search(r"<title[^>]*>(.*?)</title>", html, re.S | re.I)
            title = title_match[1].strip() if title_match else text.splitlines()[0][:160]
            project.update({"page_source": {"repo": repo, "commit": ref, "path": path,
                                             "url": f"https://github.com/buffedlizard55-lab/{repo}/blob/{ref}/{path}",
                                             "sha256": hashlib.sha256(response.stdout).hexdigest()},
                            "title": title, "opening_excerpt": text[:2200],
                            "method_terms_observed": [word for word in ("Euler", "gradient", "curvature", "dotted", "prun", "worm", "SGMC", "holdout", "HOLD", "PINN", "transfer") if re.search(r"\b" + re.escape(word) + (r"\w*" if word in ("prun", "worm", "gradient") else "") + r"\b", text, flags=re.I)],
                            "tif_links_in_page": [u for u in parser.links if ".tif" in u.lower()],
                            "status": "pinned page source read; scores below remain owner reports"})
            break
        else:
            project["status"] = "no readable page source; no invented result"
        project["source_attempts"] = attempts
        project["audited_unique_tiff_blobs"] = sum(any(a["repo"] == repo for a in e["artifacts"]) for e in inventory["unique_rasters"])
        print(repo, len(project["submissions"]), project["status"], flush=True)
        return project
    with ThreadPoolExecutor(max_workers=5) as pool:
        reviewed = list(pool.map(review, projects))
    report = {"reviewed_utc": datetime.now(timezone.utc).isoformat(),
              "scope": "Every project/score entry listed in the user brief, including empty 40–44 placeholders; GitHub source check, not live page uptime",
              "score_policy": "User reports are preserved exactly. An account-level public leaderboard row does not verify the TIFF attribution. Null means unknown, never zero.",
              "source_prompt_sha256": hashlib.sha256((ROOT / "docs/user-prompt-20261006.md").read_bytes()).hexdigest(),
              "projects": reviewed, "project_count": len(reviewed),
              "reported_score_count": sum(s["score"] is not None for p in reviewed for s in p["submissions"])}
    (ROOT / "docs/data/project-site-review-20261006.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
