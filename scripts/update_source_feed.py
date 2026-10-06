#!/usr/bin/env python3
"""Refresh allowed official-source context. NEVER polls DrivenData.

Runs on the Pages build/schedule, with TLS verification and bounded requests.
A failed update preserves last-good context while clearly marking it stale.
This feed cannot promote a candidate or alter the scientific raster.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode, urlsplit
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
PUBLIC_FEED_URL = "https://buffedlizard55-lab.github.io/GEMSDOE40/docs/data/source-feed.json"
ALLOWED_HOSTS = {"earthquake.usgs.gov", "www.usgs.gov", "www.sciencebase.gov", "gdr.openei.org"}
SOURCES = [
    ("GeoDAWN USGS metadata", "https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and"),
    ("INGENIOUS official collection", "https://gdr.openei.org/submissions/1391"),
    ("Regional geophysical metadata", "https://www.sciencebase.gov/catalog/item/628d4fabd34ef70cdba3c4a4?format=json"),
]


def earthquake_url(now: datetime) -> str:
    params = {"format": "geojson", "starttime": (now-timedelta(days=30)).date().isoformat(),
              "endtime": now.strftime("%Y-%m-%dT%H:%M:%S"),
              "minlatitude": 37.33119, "maxlatitude": 40.72788,
              "minlongitude": -120.03717, "maxlongitude": -116.14092,
              "minmagnitude": 1, "orderby": "time", "limit": 20}
    return "https://earthquake.usgs.gov/fdsnws/event/1/query?" + urlencode(params)


def get(url: str) -> tuple[bytes, int]:
    if urlsplit(url).hostname not in ALLOWED_HOSTS:
        raise ValueError("feed is restricted to the explicitly allowed official sources")
    req = Request(url, headers={"User-Agent": "GEMSDOE40-official-context/1.0 (+https://github.com/buffedlizard55-lab/GEMSDOE40)"})
    with urlopen(req, timeout=20) as response:
        if urlsplit(response.url).hostname not in ALLOWED_HOSTS:
            raise ValueError("unexpected redirect outside official-source allowlist")
        data = response.read(2_000_001)
        if len(data) > 2_000_000:
            raise ValueError("source response exceeded 2 MB safety limit")
        return data, response.status


def parse_events(data: bytes) -> list[dict]:
    document = json.loads(data)
    if document.get("type") != "FeatureCollection" or not isinstance(document.get("features"), list):
        raise ValueError("unexpected USGS GeoJSON schema")
    events = []
    for feature in document["features"]:
        props = feature["properties"]
        url = props.get("url", "")
        if urlsplit(url).hostname != "earthquake.usgs.gov":
            raise ValueError("unexpected event-link host")
        events.append({"id": feature["id"], "magnitude": props.get("mag"), "place": props.get("place"),
                       "time_ms": props.get("time"), "url": url})
    return events


def update(output: Path, restore_live: bool = False) -> dict:
    now = datetime.now(timezone.utc)
    old = json.loads(output.read_text()) if output.exists() else {}
    restore_status = "not requested"
    if restore_live:
        # A fresh Actions checkout would otherwise lose the last successful
        # generated feed. Read only this project's OWN public snapshot.
        try:
            request = Request(PUBLIC_FEED_URL, headers={"User-Agent": "GEMSDOE40-feed-preservation/1.0"})
            with urlopen(request, timeout=20) as response:
                previous = json.loads(response.read(2_000_000))
            if previous.get("last_success_utc") and (not old.get("last_success_utc") or previous["last_success_utc"] > old["last_success_utc"]):
                old = previous
            restore_status = "own deployed snapshot read"
        except Exception as exc:
            restore_status = "unavailable: " + str(exc)[:300]
    feed = {"attempted_utc": now.isoformat(), "last_success_utc": old.get("last_success_utc"),
            "status": "updating", "scope": "USGS/GDR context only, not a fault label or resource discovery feed",
            "leaderboard_monitoring": "disabled; no permission to automate DrivenData monitoring",
            "events": old.get("events", []), "source_status": [], "errors": [],
            "own_previous_snapshot_restore": restore_status}
    query = earthquake_url(now)
    feed["earthquake_query"] = query
    try:
        payload, status = get(query)
        feed["events"] = parse_events(payload)
        feed["events_sha256"] = hashlib.sha256(payload).hexdigest()
        feed["last_success_utc"] = now.isoformat()
        feed["status"] = "ok"
    except Exception as exc:
        feed["status"] = "stale" if feed["last_success_utc"] else "unavailable"
        feed["errors"].append({"source": "USGS ComCat", "error": str(exc)[:500]})
    for name, url in SOURCES:
        try:
            payload, status = get(url)
            feed["source_status"].append({"name": name, "url": url, "http_status": status,
                                         "bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()})
        except Exception as exc:
            feed["source_status"].append({"name": name, "url": url, "http_status": None, "error": str(exc)[:500]})
    feed["source_probe_failures"] = sum(s["http_status"] is None for s in feed["source_status"])
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".tmp")
    temporary.write_text(json.dumps(feed, indent=2, allow_nan=False) + "\n")
    temporary.replace(output)
    return feed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "docs/data/source-feed.json")
    parser.add_argument("--restore-live", action="store_true", help="preserve last-good context from this project's own deployed snapshot on an ephemeral runner")
    args = parser.parse_args()
    result = update(args.output, restore_live=args.restore_live)
    print(f"Official context feed: {result['status']}; events={len(result['events'])}; last_success={result['last_success_utc']}")
    # A blocked network must not hide the verified artifact by failing Pages.


if __name__ == "__main__":
    main()
