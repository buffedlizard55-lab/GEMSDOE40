#!/usr/bin/env python3
"""Verify real public Pages routing and download bytes after deployment.

No simulated success and no competition endpoints. Intended for the Pages
Actions runner (sandbox direct HTTPS may be unavailable). Bounded retries allow
CDN propagation; failure is a failing workflow with the actual error.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
BASE = "https://buffedlizard55-lab.github.io/GEMSDOE40/"


def get(relative):
    request = Request(BASE + relative, headers={"User-Agent": "GEMSDOE40-deployment-integrity/1.0", "Cache-Control": "no-cache"})
    with urlopen(request, timeout=30) as response:
        return response.read(8_000_000), response.status, response.headers.get("Content-Type")


def main():
    current = json.loads((ROOT / "docs/data/current-candidate.json").read_text())
    expected_html = hashlib.sha256((ROOT / "docs/index.html").read_bytes()).hexdigest()
    last_error = None
    for attempt, delay in enumerate((0, 5, 10, 20, 30, 45), 1):
        if delay: time.sleep(delay)
        try:
            rows = []
            for path in ("docs/index.html", "docs/executive-summary.html", "executive-summary.html"):
                body, status, mime = get(path)
                assert status == 200 and current["filename"].encode() in body, f"stale/missing current download on {path}"
                if path == "docs/index.html": assert hashlib.sha256(body).hexdigest() == expected_html, "old index HTML still at CDN"
                rows.append({"url": BASE+path, "status": status, "content_type": mime})
            root, status, mime = get("")
            assert status == 200 and b"docs/index.html" in root, "root redirect route missing"
            for prefix in ("docs/downloads/", "downloads/"):
                body, status, mime = get(prefix + current["filename"])
                digest = hashlib.sha256(body).hexdigest()
                assert status == 200 and digest == current["sha256"], f"served TIFF SHA-256 mismatch at {prefix}"
                rows.append({"url": BASE+prefix+current["filename"], "status": status, "bytes": len(body), "content_type": mime, "sha256": digest})
            feed_body, status, _ = get("docs/data/source-feed.json")
            feed = json.loads(feed_body)
            assert status == 200 and feed["status"] in ("ok", "stale", "unavailable")
            report = {"verified_utc": datetime.now(timezone.utc).isoformat(), "pass": True,
                      "attempt": attempt, "checks": rows, "source_feed_status": feed["status"],
                      "source_feed_last_success": feed.get("last_success_utc")}
            out = ROOT / "work/deployment-verification.json"; out.parent.mkdir(exist_ok=True)
            out.write_text(json.dumps(report, indent=2)+"\n")
            print(json.dumps(report, indent=2))
            return
        except Exception as exc:
            last_error = exc
            print(f"Public verification attempt {attempt}/6: {type(exc).__name__}: {exc}", flush=True)
    raise RuntimeError(f"Cannot claim verified public deployment: {last_error}")


if __name__ == "__main__":
    main()
