#!/usr/bin/env python3
"""Exercise actual desktop/mobile pages, download bytes, search and copy controls.

Optional QA dependency: playwright and Chromium. Uses an existing --url or a
short-lived local server for the already staged _site directory. No competition
network access or upload occurs.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
from functools import partial
import hashlib
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
import os
from pathlib import Path
import threading

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
# The five pages the current generator actually builds (methodology / leaderboard pages were
# retired with the H4-era site; their receipts remain, their routes do not).
PAGES = ("index.html", "executive-summary.html", "evidence.html", "hypotheses.html", "sources.html")


class QuietHandler(SimpleHTTPRequestHandler):
    def log_message(self, *args): pass


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", help="existing local site URL, with its current-page base path")
    parser.add_argument("--output", type=Path, default=ROOT / "work/browser-checks.json")
    args = parser.parse_args()
    server = None
    if args.url:
        base = args.url.rstrip("/") + "/"
    else:
        server = ThreadingHTTPServer(("0.0.0.0", 0), partial(QuietHandler, directory=str(ROOT/"_site")))
        threading.Thread(target=server.serve_forever, daemon=True).start()
        base = f"http://127.0.0.1:{server.server_port}/docs/"
    release = json.loads((ROOT / "docs/data/current-candidate.json").read_text())
    checks, errors = [], []
    snapshots = ROOT / "work/browser-screenshots"; snapshots.mkdir(parents=True, exist_ok=True)
    try:
        with sync_playwright() as p:
            executable = os.getenv("CHROMIUM_EXECUTABLE_PATH")
            browser = p.chromium.launch(headless=True, executable_path=executable,
                                       args=["--no-sandbox", "--disable-dev-shm-usage", "--single-process", "--no-zygote"])
            context = browser.new_context(accept_downloads=True, permissions=["clipboard-read", "clipboard-write"])
            page = context.new_page()
            page.on("pageerror", lambda error: errors.append(str(error)))
            for width, height in ((1440, 1000), (390, 844), (320, 800)):
                page.set_viewport_size({"width": width, "height": height})
                for name in PAGES:
                    response = page.goto(base + name, wait_until="networkidle")
                    assert response and response.status == 200, (width, name, "HTTP")
                    assert not page.evaluate("document.documentElement.scrollWidth > window.innerWidth"), (width, name, "body overflow")
                    assert page.locator("h1").count() == 1, (name, "one main title required")
                    assert page.locator("a.skip").get_attribute("href") == "#main"
                    if name in ("index.html", "executive-summary.html"):
                        link = page.locator(f'a[download][href$="{release["filename"]}"]').first
                        assert link.is_visible()
                        if name == "index.html":
                            box = link.bounding_box()
                            assert box and box["y"] < height, (width, "primary download begins below first viewport")
                    checks.append({"page": name, "width": width, "http": 200, "no_horizontal_body_overflow": True})
                    if name == "index.html": page.screenshot(path=str(snapshots/f"overview-{width}.png"), full_page=True)
            page.set_viewport_size({"width": 1440, "height": 1000})
            page.goto(base + "executive-summary.html", wait_until="networkidle")
            for name in ("filename", "submission-name", "submission-note", "file-sha"):
                page.locator(f'[data-copy="{name}"]').click()
                page.wait_for_function("id => document.querySelector(`[data-copy=\"${id}\"]`).textContent === 'Copied'", arg=name)
                actual = page.evaluate("navigator.clipboard.readText()")
                assert actual == page.locator(f"#{name}").text_content().strip(), name
            with page.expect_download() as pending:
                page.locator(f'a[download][href$="{release["filename"]}"]').click()
            download = pending.value
            assert download.suggested_filename == release["filename"]
            actual_sha = hashlib.sha256(Path(download.path()).read_bytes()).hexdigest()
            assert actual_sha == release["sha256"], "browser downloaded different bytes"
            page.goto(base + "leaderboard.html", wait_until="networkidle")
            page.locator("#project-search").fill("GEMSDOE32")
            assert page.locator("[data-project-row]:visible").count() == 1
            page.locator("#project-search").fill("")
            assert page.locator("[data-project-row]:visible").count() == 44
            page.locator("#project-search").fill("no-project-with-this-name")
            assert page.locator("[data-project-row]:visible").count() == 0
            assert not errors, errors
            browser.close()
    finally:
        if server: server.shutdown(); server.server_close()
    receipt = {"verified_utc": datetime.now(timezone.utc).isoformat(), "pass": True,
               "pages_and_viewports": checks, "download_name": release["filename"],
               "download_sha256": actual_sha, "copy_controls_verified": 4,
               "project_search_verified": True, "javascript_errors": errors,
               "scope": "Local Chromium browser on staged/static public files; live deployment checked separately"}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"BROWSER PASS: {len(checks)} page/viewport combinations; exact downloaded TIFF, 4 copy controls, project search; no JS errors")


if __name__ == "__main__":
    main()
