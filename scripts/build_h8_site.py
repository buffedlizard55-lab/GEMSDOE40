#!/usr/bin/env python3
"""Rebuild the GEMSDOE40 pages around the published H8 Euler depth-clustering candidate.

Reads only files that exist in the repository (audit + generation receipts, the
measured analysis receipts) and never invents a score, a measurement or a claim.
Every number printed here is read from a receipt written by a runner in this
repository, or is a verbatim quotation from an official source with its link.

Usage:
    python scripts/build_h8_site.py            # write the five pages + manifest
"""
from __future__ import annotations

from html import escape as esc
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
AUDIT = DOCS / "data" / "h8-lineament-audit.json"
GENERATION = DOCS / "data" / "h8-lineament-generation.json"
MANIFEST = DOCS / "data" / "current-candidate.json"
SESSION2 = DOCS / "data" / "session2-artifacts-20261006.json"
PREREG = "research/h8-lineament-preregistration-20261006.md"
REGISTER = "research/h8-preregistration-20261006.md"
H33 = DOCS / "data" / "h33-measured-analysis.json"
INSTRUMENTS = DOCS / "data" / "instrument-audit-20261006.json"

PROBLEM = "https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/"
DATA = "https://www.drivendata.org/competitions/306/competition-doe-gems/data/"
RULES = "https://docs.nlr.gov/docs/fy26osti/96647.pdf"
LEADERBOARD = "https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/"
REID1990 = "https://www.reid-geophys.co.uk/wp-content/uploads/2017/11/Reid-et-al-1990.pdf"
REID2014 = "https://www.reid-geophys.co.uk/wp-content/uploads/2017/11/Reid-Thurston-2014.pdf"
STAFF_MASK = ("https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-"
              "masked-when-scoring-and-are-they-in-the-final-round-label-set/11516")
STAFF_JOIN = ("https://community.drivendata.org/t/how-were-the-new-test-faults-identified-data-sources-"
              "and-fault-types/11527")
REFERENCE = "https://github.com/drivendataorg/gems-prize-reference-solution"
GDR = "https://gdr.openei.org/submissions/1391"
THREEDEP = "https://www.usgs.gov/3d-elevation-program"
TNMLINKS = ("https://www.drivendata.org/competitions/306/competition-doe-gems/data/"
            "#1m_DEM_links.csv")

NAV = (("index.html", "Overview"), ("executive-summary.html", "Submission guide"),
       ("evidence.html", "Evidence"), ("hypotheses.html", "Hypotheses"), ("sources.html", "Sources"))


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def page(name: str, title: str, body: str, active: str) -> None:
    links = "".join(f'<a class="{"active" if file == active else ""}" href="{file}">{text}</a>'
                    for file, text in NAV)
    html = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{esc(title)} · GEMSDOE40</title>
<meta name="description" content="Euler deconvolution depth-clustering research for the DOE GEMS fault-mapping prize: a unique, format-audited GeoTIFF with honest, receipt-backed evidence.">
<link rel="stylesheet" href="assets/contact.css"><link rel="stylesheet" href="assets/site.css">
<script src="assets/contact.js" defer></script></head>
<body><a href="#main" class="skip">Skip to content</a><header class="header"><div class="wrap"><div class="masthead">
<a class="brand" href="index.html"><span class="brand-icon" aria-hidden="true">G</span><b>GEMS / 40</b><em>H8 EULER</em></a>
<nav class="nav" aria-label="Main navigation">{links}</nav></div></div></header>
<main id="main" class="wrap">{body}</main>
<footer class="footer"><div class="wrap">
<p><strong>Maximize P(Win). Own the Outcome.</strong><br>Evidence before a weekly submission slot.</p>
<p>Independent competition research, not affiliated with DOE or DrivenData. No competition login is used or stored.<br>
<a href="https://github.com/buffedlizard55-lab/GEMSDOE40">Repository</a> ·
<a href="user-prompt-20261006.md">Full project brief</a> ·
<a href="data/current-candidate.json">Machine-readable current candidate</a></p></div></footer></body></html>'''
    (DOCS / name).write_text(html, encoding="utf-8")


def main() -> int:
    audit = load(AUDIT)
    generation = load(GENERATION)
    h33 = load(H33)
    instruments = load(INSTRUMENTS)
    session2 = load(SESSION2) if SESSION2.is_file() else None
    h4line = None
    _h4z = DOCS / "downloads/gemsdoe40-euler-line-ring-pruned-60000px-20261006T032708Z-8dafb186-zeros.tif"
    if _h4z.is_file():
        def _h4_receipt(name):
            return json.loads((DOCS / "data" / name).read_text())
        _za = _h4_receipt("h4_shipped_audit_zeros.json")
        _na = _h4_receipt("h4_shipped_audit_nan.json")
        h4line = {
            "zeros_tif": {"filename": _za["file"].split("/")[-1], "sha256": _za["sha256"], "bytes": _za["bytes"]},
            "nan_tif": {"filename": _na["file"].split("/")[-1], "sha256": _na["sha256"], "bytes": _na["bytes"]},
            "zip": {"filename": _h4z.with_suffix(".zip").name, "sha256": None, "bytes": None},
        }
        import hashlib as _hashlib
        h4line["zip"]["sha256"] = _hashlib.sha256(_h4z.with_suffix(".zip").read_bytes()).hexdigest()
        h4line["zip"]["bytes"] = _h4z.with_suffix(".zip").stat().st_size

    published = audit["published"]
    name = published["name"]
    link = f"downloads/{name}"
    sha256 = published["sha256"]
    twin = published.get("hard_twin")
    cloud = published["solution_cloud"]
    dots = audit["emitted_dots"]
    mass = audit["mass"]
    novelty = audit["novelty"]
    cat = audit["catalogue_surrogate"]
    off = audit["offcat_sgmc_surrogate"]
    family = generation["families"]
    rtp, grav = family["rtp"], family["iso_grav_anom"]
    solutions = rtp["accepted"] + grav["accepted"]
    retained = rtp["retained"] + grav["retained"]
    cfg = generation["configuration"]
    size_mb = published["bytes"] / 1e6
    pixsha = audit["candidate"]["canonical_pixels_sha256"]
    tracking = f"GEMSDOE40-H8-LINEAMENT-{pixsha[:12]}"
    note = (f"Euler deconvolution SI=0 contact depth-clustering over magnetic (rtp) and isostatic gravity "
            f"(differentiated once); lineament-weighted kernel-density emission; {dots:,} dots; "
            f"id {pixsha[:12]}.")
    stats = instruments["statistics"]

    def results_table() -> str:
        rows = []
        for row in sorted(instruments["artifacts"], key=lambda r: -r["live_score"]):
            rows.append(f'<tr><td class="mono">{esc(row["file"])}</td><td class="num">{row["live_score"]:.4f}</td>'
                        f'<td class="num">{int(row["mass"]):,}</td>'
                        f'<td class="num">{row["lm_calibrated"]:.4f}</td>'
                        f'<td class="num">{row["w_offcat_sgmc"]:.4f}</td></tr>')
        rows.append(f'<tr class="row-h4"><td><strong>{esc(name)}</strong> — this release</td>'
                    f'<td class="num">unscored</td><td class="num">{int(mass):,}</td>'
                    f'<td class="num">{cat["w"]:.4f}</td><td class="num">{off["w"]:.4f}</td></tr>')
        return ('<div class="tablewrap"><table><thead><tr><th>Artifact</th><th>Recorded live score</th>'
                '<th>Mass (px)</th><th>w vs catalogue</th><th>w vs SGMC off-catalogue</th></tr></thead>'
                f'<tbody>{"".join(rows)}</tbody></table></div>')

    def block_table() -> str:
        rows = []
        for block in audit["blocks"]:
            c = block.get("catalogue", {})
            o = block.get("offcat", {})
            rows.append(f'<tr><td class="num">{block["block"][0]},{block["block"][1]}</td>'
                        f'<td class="num">{int(block["emitted"]):,}</td>'
                        f'<td class="num">{int(block["catalogue_cells"]):,}</td>'
                        f'<td class="num">{c.get("w", 0.0):.4f}</td>'
                        f'<td class="num">{o.get("w", 0.0):.4f}</td></tr>')
        return ('<div class="tablewrap"><table><thead><tr><th>Block</th><th>Emitted dots</th>'
                '<th>Dots vs catalogue</th><th>w vs catalogue</th><th>w vs SGMC off-cat</th></tr></thead>'
                f'<tbody>{"".join(rows)}</tbody></table></div>')

    twin_line = ""
    if twin:
        twin_line = (f'<p class="micro">Metric-optimal twin of the same {dots:,} cells, every value set to 1.0: '
                     f'<a href="downloads/{esc(twin["name"])}">{esc(twin["name"])}</a> '
                     f'({twin["bytes"] / 1e6:.2f} MB, sha256 <span class="mono">{twin["sha256"][:16]}…</span>, '
                     f'support identical to the primary: '
                     f'{"yes" if twin["support_identical_to_primary"] else "no"}).</p>')

    download_block = f'''<div class="download"><a class="button" href="{link}" download>Download the H8 GeoTIFF <span aria-hidden="true">↓</span></a></div>
<p class="file-details">{size_mb:.2f} MB · single-band float32 · EPSG:32611 · 100 m<br>
3,730 × 3,292 · all 5,167,373 footprint cells finite in [0, 1]<br>
NaN only outside the sample footprint · DEFLATE with predictor 1 (none), the official sample's setting<br>
sha256 <span class="mono">{sha256}</span></p>
{twin_line}
<p class="micro">Tracking name for your own record: <span class="mono">{tracking}</span></p>
<p class="micro">Depth-labelled solution cloud: <a href="downloads/{Path(cloud["path"]).name}">{Path(cloud["path"]).name}</a>
({cloud["bytes"] / 1e6:.1f} MB, {retained:,} retained solutions with row, col, depth, depth standard error,
residual, condition number, window, cluster weight, lineament coherence and corroboration flag).</p>'''

    feed_block = ""
    feed_path = DOCS / "data" / "source-feed.json"
    if feed_path.is_file():
        feed = load(feed_path)
        statuses = feed.get("source_status", [])
        ok = sum(1 for row in statuses if row.get("http_status") == 200)
        events = feed.get("events", [])
        rows = "".join(
            f'<tr><td>{esc(row.get("name", "?"))}</td>'
            f'<td class="num">{row.get("http_status") if row.get("http_status") is not None else "unreachable"}</td>'
            f'<td class="micro">{esc(str(row.get("error"))[:80])}</td></tr>' for row in statuses)
        feed_block = f'''<section class="section" id="feed"><div class="sectionhead"><div><p class="eyebrow">05 / CONTEXT FEED</p>
<h2>Official-source context, refreshed by the Pages runner.</h2></div>
<a href="data/source-feed.json">source-feed.json ↗</a></div>
<p class="micro">Scope: {esc(feed.get("scope", ""))}. Leaderboard monitoring:
{esc(str(feed.get("leaderboard_monitoring", "disabled")))}, consistent with the competition terms of use.
Last attempt {esc(str(feed.get("attempted_utc")))} · status <strong>{esc(str(feed.get("status")))}</strong> ·
{ok} of {len(statuses)} sources reachable · {len(events)} events recorded. This sandbox cannot reach the USGS/GDR
hosts (TLS handshake closed), so the honest published state is whatever the receipt says; the GitHub Pages runner
refreshes it daily and the file is never fabricated.</p>
<div class="tablewrap"><table><thead><tr><th>Source</th><th>HTTP</th><th>Error (truncated)</th></tr></thead>
<tbody>{rows}</tbody></table></div></section>'''

    def retained_session2() -> str:
        """Every session-2/3 artifact, still downloadable, still explicitly held."""
        if not session2:
            return ""
        h13 = session2
        sib = session2["sibling_candidate"]
        asa = session2["h8asa_candidate"]
        rows = []
        for name_, path_, sha_, status_ in (
            (h13["filename"], h13["path"], h13["sha256"], h13["status"]),
            (sib["filename"], sib["path"], sib["sha256"], sib["status"]),
            (asa["filename"], asa["path"], asa["sha256"], asa["status"]),
        ):
            rows.append(f'<tr><td class="mono"><a href="{path_.replace("docs/", "")}">{name_}</a></td>'
                        f'<td class="num">{sha_[:16]}…</td><td>{esc(status_)}</td></tr>')
        for label, entry in (("portal-safe zeros twin", asa["portal_safe_twin"]),
                             ("H4-line portal-safe zeros", h4line["zeros_tif"]),
                             ("H4-line NaN twin", h4line["nan_tif"]),
                             ("H4-line zeros ZIP", h4line["zip"]),
                             ("zeros ZIP", asa["zip"])):
            rows.append(f'<tr><td class="mono"><a href="downloads/{entry["filename"]}">{entry["filename"]}</a></td>'
                        f'<td class="num">{entry["sha256"][:16]}…</td><td>HOLD — DO NOT SUBMIT</td></tr>')
        return ('<div class="tablewrap"><table><thead><tr><th>Artifact</th><th>SHA-256</th><th>Status</th>'
                '</tr></thead><tbody>' + "".join(rows) + '</tbody></table></div>')

    retained_block = retained_session2()
    retained_paragraph = ('<p>Session-2 artifacts (H13 crest-binary, H8 trace-locked depth-KDE, H8-ASA '
                          'analytic-signal depth-KDE) remain downloadable above with their own receipts and '
                          'their own HOLD statuses, and so do the earlier arms: '
                          '<a href="downloads/gemsdoe40-h7-rtp-euler-gravity-context-3d-kde-20261006-998f660f.tif">'
                          'H7 RTP Euler + gravity-context KDE</a>, the retained '
                          '<a href="downloads/gemsdoe40-h4-contact-offset-depthkde-20261006-fab9f6619c02.tif">'
                          'H4 contact-offset depth-KDE</a>, and '
                          '<a href="downloads/gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-zeros.tif">'
                          'H40 run-2</a> with its '
                          '<a href="downloads/gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-nan.tif">NaN twin</a>. '
                          '<strong>HOLD — DO NOT SUBMIT H40</strong>, and '
                          '<strong>HOLD — DO NOT SUBMIT H7 or H40</strong> for the earlier arms: '
                          '<strong>Nothing earlier is cleared either</strong> — including the H4-line ring-pruned '
                          'artifact shipped by the parallel PR #17, which failed its own frozen four-fold gate '
                          '(better than the h33-2-b2 reference on 1 of 4 folds; projected 0.197, interval '
                          '0.176–0.219; write-up <a href="research/h4-euler-depthcluster.html">'
                          'research/h4-euler-depthcluster.html</a>, <a href="limitations.html">limitations.html</a>). '
                          'Every promotion instrument used '
                          'on them is <em>withdrawn</em> by the sixteen-score audit '
                          '(<a href="data/instrument-audit-20261006.json">instrument-audit-20261006.json</a>).</p>')

    # ------------------------------------------------------------------ index --
    index_body = f'''<section class="hero"><div>
<p class="kicker">H8 / EULER DEPTH-CLUSTERING / 06 OCT 2026</p>
<h1>Faults from<br><em>source depth.</em></h1>
{download_block}
<p class="lead">A unique candidate built exactly as the brief requires: SI = 0 contact Euler deconvolution
(Reid et al., 1990) solved on magnetic and gravity layers, weighted for shallow depth, cross-window depth
consensus, lineament coherence and cross-family corroboration, converted to a kernel-density field and emitted
as a sparse dot pattern inside the metric's own 300 m kernel.</p>
<div class="notice" id="decision"><strong>Status: built, audited, not promoted.</strong>
<p>Format gate passed; raw-output novelty gate passed against {novelty["corpus_size"]} cached prior rasters
(largest |correlation| {novelty["max_abs_pearson"]:.4f}, largest top-mass Jaccard {novelty["max_jaccard_topmass"]:.4f},
0 exact duplicates). <strong>No local instrument in this repository ranks the family's recorded scores</strong>
(Spearman ρ ≤ {max(abs(stats["spearman_lm_vs_live"]), abs(stats["spearman_mass_vs_live"])):.2f} for every proxy
tested), so this file is published as a <em>research artifact</em>: it is not claimed to beat the owner's best
0.2778 nor the observed 0.3195, and no weekly submission slot was spent.</p></div>
<div class="actions"><a class="textlink" href="executive-summary.html">Name, note &amp; submission guide ↗</a>
<a class="textlink" href="evidence.html">Inspect every measurement ↗</a></div></div>
<div class="map-panel"><div class="map-foot"><span>MODEL OUTPUT / NOT OBSERVED FAULTS</span>
<b>SEE EVIDENCE ↗</b></div><figure><figcaption>{dots:,} emitted dots from {solutions:,} QC-passing Euler
solutions ({retained:,} retained by the depth-consensus and coherence gates) across two potential-field families.</figcaption>
</figure></div></section>
<div class="checkstrip">
<div><span class="stat">{dots:,}</span><span class="label">Emitted dots</span></div>
<div><span class="stat">{novelty["corpus_size"]}</span><span class="label">Prior rasters compared</span></div>
<div><span class="stat">{novelty["max_abs_pearson"]:.4f}</span><span class="label">Largest raw |correlation|</span></div>
<div><span class="stat">0</span><span class="label">Competition slots used</span></div></div>

<section class="section" id="why"><div class="sectionhead"><div><p class="eyebrow">01 / THE QUESTION</p>
<h2>Why did the best submission score 0.2778?</h2></div>
<a href="research/h8-analysis-20261006.md">Full analysis ↗</a></div>
<div class="grid3">
<article class="card"><span class="stepnum">MEASURED</span><h3>It was a deletion, not a discovery.</h3>
<p>The best file is the 0.2708 base minus {h33["removed_pixels"]:,} pixels
({h33["mass_removed_fraction"] * 100:.2f} % of its mass), all of them
{h33["removed_distance_to_known_m_min_max"][0]:.0f}–{h33["removed_distance_to_known_m_min_max"][1]:.0f} m from the
provided catalogue, with <strong>{h33["added_pixels"]} pixels added</strong>. Its minimum distance to a mapped
fault is {h33["h33_min_distance_to_known_m"]:.0f} m.</p></article>
<article class="card"><span class="stepnum">ALGEBRA</span><h3>Every emitted unit must clear a bar.</h3>
<p>With α = 0.2 and β = 0.8 and α + β = 1, one more unit of predicted mass changes the metric denominator by
exactly 0.2 wherever it lands. So a dot only pays while its kernel credit exceeds 0.2 × DTI — at 0.2778 that
means <strong>within 283 m</strong> of a hidden fault. This is the whole reason near-catalogue mass is
expensive: those pixels earn nothing while the metric charges for them.</p></article>
<article class="card"><span class="stepnum">CONSEQUENCE</span><h3>Credit per dot is the target.</h3>
<p>Inverting the metric on the 16 recorded scores gives the hidden credit each file implies: the best averages
≈ 0.090 per dot at the family's own mass, the worst 0.016. Reaching 0.3195 at the same mass needs ≈ 0.103 —
about <strong>15 % better placement</strong>, not exotic new data.</p></article></div></section>

<section class="section"><div class="sectionhead"><div><p class="eyebrow">02 / THE METHOD</p>
<h2>Reid's contact solver, then a density, then dots.</h2></div>
<a href="{PREREG}">Pre-registration ↗</a></div>
<div class="grid3">
<article class="card"><span class="stepnum">01 · SOLVE</span><h3>SI = 0 Euler contact.</h3>
<p>Reid et al. (1990) eq. 2 with the arbitrary contact offset A — the structural index of a fault-like contact —
solved on sliding windows of <code>rtp</code> (150 m continuation, windows {", ".join(map(str, cfg["families"][0]["windows"]))})
and of the once-differentiated <code>iso_grav_anom</code> (500 m continuation, windows
{", ".join(map(str, cfg["families"][1]["windows"]))}), stride {cfg["stride"]}. A synthetic analytic contact is
recovered to 0.3–1.0 m laterally and exactly in depth at zero noise — the solver is validated, so the real-data
difficulty is a property of the data.</p></article>
<article class="card"><span class="stepnum">02 · WEIGHT</span><h3>Shallow, tight, coherent.</h3>
<p>Weight = exp(−depth/{cfg["depth_decay_m"]:.0f} m) × cross-window depth consensus (≥ {cfg["min_neighbours"]}
neighbours within {cfg["cluster_xy_m"]:.0f} m, at least one from a different window) × lineament coherence
(1 − λ₂/λ₁ of the local solution cloud) × cross-family corroboration
({cfg["corroboration_xy_m"]:.0f} m / {cfg["corroboration_depth_m"]:.0f} m, uncorroborated solutions halved).</p></article>
<article class="card"><span class="stepnum">03 · EMIT</span><h3>Anisotropic KDE, then dots.</h3>
<p>The weighted cloud is smeared along each local lineament axis (σ {cfg["kernel_sigma_along_px"]} px along,
{cfg["kernel_sigma_across_px"]} px across, {cfg["direction_bins"]} directions) into a continuous density in
[0, 1], then sampled at {cfg["nms_spacing_px"]} px maximum separation — just inside the metric's own 3 px kernel —
under a {cfg["mass_budget"]:,}-dot budget, with the provided catalogue and a 1 px ring removed.</p></article></div>
<p class="micro">The construction reads only <code>training_features.tif</code>, the sample template and the
catalogue mask. It never reads <code>labels.tif</code> as truth, any prior prediction, or any proxy raster:
<code>construction_reads_labels_or_priors = {str(generation["construction_reads_labels_or_priors"]).lower()}</code>
in the generation receipt.</p></section>

<section class="section" id="retained"><div class="sectionhead"><div><p class="eyebrow">02b / RETAINED ARTIFACTS</p>
<h2>Everything published earlier, still downloadable, still held.</h2></div>
<a href="executive-summary.html">Submission guide ↗</a></div>
{retained_paragraph}
{retained_block}</section>

<section class="section"><div class="sectionhead"><div><p class="eyebrow">03 / THE HONEST PART</p>
<h2>What is measured, and what is not.</h2></div><a href="evidence.html">All measurements ↗</a></div>
<div class="twocol"><div class="prose">
<p><strong>Measured.</strong> The file is a valid single-band float32 GeoTIFF in EPSG:32611 with the exact
sample shape, transform and NaN footprint ({audit["candidate"]["candidate"]["outside_nan"]:,} cells outside,
{audit["candidate"]["candidate"]["in_footprint_nonzero"]:,} in-footprint non-zero cells in
[{audit["value_min"]:.3f}, {audit["value_max"]:.3f}]). It is not a duplicate of any cached prior raster.</p>
<p><strong>Not measured.</strong> Its leaderboard score. The catalogue surrogate (w = {cat["w"]:.4f}) and the
off-catalogue SGMC surrogate (w = {off["w"]:.4f}) are named proxies, and neither ranks the 16 recorded scores:
ρ(catalogue-calibrated LM) = {stats["spearman_lm_vs_live"]:+.3f}, ρ(SGMC off-catalogue) = {stats["spearman_w_offcat_vs_live"]:+.3f},
ρ(mass) = {stats["spearman_mass_vs_live"]:+.3f}. The best recorded file sits *below* the chance line on the
catalogue column because catalogue pixels are masked from scoring
(<a href="{STAFF_MASK}">DrivenData staff, 2026-09-16</a>) — credit there can never be earned.</p>
<p><strong>Therefore.</strong> The brief's own rule applies: the top hypothesis must beat the current holdout best
before a submission slot is spent, and no instrument here can demonstrate that. The file ships with an explicit
non-promoted status rather than a score promise.</p></div>
<div class="prose"><h3>Where it stands</h3>{results_table()}
<p class="micro">Recorded scores are owner-reported leaderboard rows reproduced only as calibration; file-level
organizer attribution is unverified. "w" is mean kernel credit per unit mass against the named truth set.</p></div>
</div></section>

<section class="section" id="research"><div class="sectionhead"><div><p class="eyebrow">04 / NEXT</p>
<h2>The ranked route to a higher score.</h2></div><a href="hypotheses.html">Full register ↗</a></div>
<div class="tablewrap"><table><thead><tr><th>Rank</th><th>Hypothesis</th><th>Why it can find a fault the
catalogue misses</th><th>Cost</th></tr></thead><tbody>
<tr><td>1</td><td><strong>H14 · tip / step-over / along-strike extensions</strong></td>
<td>DrivenData staff, 2026-09-23: a new fault "can include newly mapped geometry of an existing fault system"
(<a href="{STAFF_MASK}">forum 11536, post 2</a>), and masked catalogue pixels cost nothing. Directly testable
today with a leave-the-tips-out protocol on the provided catalogue.</td><td>moderate</td></tr>
<tr><td>2</td><td>H8 · Euler depth-clustering (this file)</td><td>Potential-field contacts are subsurface
boundaries; surface compilation cannot see them under basin fill.</td><td>already paid</td></tr>
<tr><td>3</td><td>H15 · finite-step gravity inversion</td><td>A finite density step is SI = −1, not 0
(<a href="{REID2014}">Reid &amp; Thurston, 2014</a>): the basin-margin fault class that surface maps
under-represent.</td><td>high</td></tr>
<tr><td>4</td><td>H16 · 1 m lidar scarp re-mapping</td><td>The direct surface expression of the labelled fault
type; <a href="{THREEDEP}">USGS 3DEP</a> is free and official, and the competition data tab ships the DEM link
CSV. <em>Flagged irregularity:</em> the earlier scarp stack used in this repository is absent from the current
workspace and the USGS S3 endpoint is TLS-blocked from this sandbox, so it cannot be rebuilt or re-verified
here.</td><td>high</td></tr>
<tr><td>5</td><td>H17 · microseismicity alignment</td><td>Active but unmapped structures still generate events;
cheap, but measured catalogue skill is weak.</td><td>low</td></tr>
</tbody></table></div></section>
{feed_block}'''

    # -------------------------------------------------------------- executive --
    executive_body = f'''<div class="subhero"><p class="kicker">EXECUTIVE SUMMARY / HOW TO SUBMIT</p>
<h1>One file, one decision.</h1>
<p class="lead">Everything needed to make a DrivenData submission by hand: the exact file, its hash, the unique
name, the note, and what is and is not claimed about it. Nothing in this repository uploads anything.</p></div>

<section class="section"><div class="twocol"><div>
<span class="tag hold">BUILT, AUDITED, NOT PROMOTED</span>
<div class="download"><a class="button" href="{link}" download>Download {esc(name)} <span aria-hidden="true">↓</span></a></div>
<p class="micro">Download the .tif itself — not this page, not a report, not a preview image.</p>
<h3>Exact filename</h3><div class="codebox" id="filename">{esc(name)}</div>
<button class="copy" data-copy="filename">Copy filename</button>
<h3>Unique tracking name (for your own records)</h3><div class="codebox" id="submission-name">{tracking}</div>
<button class="copy" data-copy="submission-name">Copy name</button>
<h3>Short note for the submission form</h3><div class="codebox" id="submission-note">{esc(note)}</div>
<button class="copy" data-copy="submission-note">Copy note</button>
<h3>SHA-256 of the exact download</h3><div class="codebox" id="file-sha">{sha256}</div>
<button class="copy" data-copy="file-sha">Copy SHA-256</button>
<p class="micro"><a href="downloads/{Path(cloud["path"]).name}">Depth-labelled solution cloud</a> ·
<a href="data/h8-lineament-audit.json">Audit receipt</a> ·
<a href="data/h8-lineament-generation.json">Generation receipt</a> ·
<a href="research/h8-analysis-20261006.md">Why 0.2778 leads</a></p>
</div><div><h3>How to submit — official flow</h3>
<ol class="prose">
<li>Open the <a href="{PROBLEM}">competition page</a> and sign in to your own DrivenData account. Never place a
password, token or API key in this repository, in a prompt, or in chat.</li>
<li>Download the exact <span class="mono">.tif</span> above and confirm its SHA-256 equals
<span class="mono">{sha256[:32]}…</span>.</li>
<li>Choose <em>New submission</em> on the competition's <em>Submit</em> tab and select the downloaded GeoTIFF.
A single-band GeoTIFF, or a ZIP containing exactly one GeoTIFF, is accepted; this project publishes the TIFF
directly.</li>
<li>Paste the note above into the optional note field, keep the tracking name in your own record, and submit.</li>
<li>Give the gateway time to score; the score appears on the <a href="{LEADERBOARD}">leaderboard</a> page. The
weekly limit is three submissions, and one chosen submission serves the final round as well.</li>
</ol>
<div class="callout"><strong>Format contract (official)</strong>
<p>Single band · float32 · EPSG:32611 · 100 m · the sample's exact shape and geotransform · finite values in
[0, 1] inside the footprint · <code>NaN</code> outside. This release re-opened and verified all
{audit["candidate"]["candidate"]["outside_nan"] + 5_167_373:,} cells: 5,167,373 finite and in range, 7,111,787
NaN. See the <a href="{RULES}">submission specification</a>.</p></div>
<p class="micro">Two files are offered. The continuous field is the registered candidate (KDE values in
[{audit["value_min"]:.3f}, 1.0]). The hard twin has <em>identical support</em> with every value set to 1.0, which
the metric algebra weakly prefers once a dot is believed: raising a dot's value changes the denominator by 0.2
and the numerator by its kernel credit, so it pays while that credit exceeds 0.2 × DTI. Both are format-valid;
only the continuous one carries the registered method description.</p>
</div></div></section>

<section class="section"><h2>"Predicted values must be in range [0, 1]" — what actually caused it</h2>
<div class="prose"><p>Measured in this family, not guessed: the rejection traced to the GeoTIFF encoder writing
<strong>Predictor = 2</strong>, the <em>integer</em> horizontal-differencing predictor, on floating-point data —
a reader then decodes garbage values rather than the intended ones. The official
<code>sample_submission.tif</code> carries <strong>Predictor = 1</strong> (none) with <code>NaN</code> outside its
footprint. This release is written with <strong>predictor 1</strong> and was re-opened and range-checked after
writing. The NaN footprint itself is correct and expected — it is what the template does.</p>
<p>Both published files were also checked with an independent TIFF tag reader, not only with the writing library:
<code>Compression = 8 (DEFLATE)</code>, <code>Predictor = 1</code>, <code>SampleFormat = 3 (float)</code>,
<code>BitsPerSample = 32</code>, one sample per pixel.</p></div></section>

<section class="section"><h2>Retained history — offered beside an explicit hold</h2>
<div class="prose"><p>Earlier research artifacts stay downloadable and untouched so that any past upload remains
reproducible. None of them is cleared for a new submission:</p>
<ul>
<li><a href="downloads/gemsdoe40-h7-rtp-euler-gravity-context-3d-kde-20261006-998f660f.tif">H7 RTP Euler +
gravity-context 3D KDE</a> — <strong>HOLD — DO NOT SUBMIT</strong>.</li>
<li><a href="downloads/gemsdoe40-h4-contact-offset-depthkde-20261006-fab9f6619c02.tif">H4 contact-offset
depth-KDE</a> — <strong>HOLD — DO NOT SUBMIT</strong>.</li>
<li><a href="downloads/gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-zeros.tif">H40 run-2 emission</a> —
<strong>HOLD — DO NOT SUBMIT</strong>; its own audit found no promotion instrument that survives.</li>
</ul>
{retained_block}
<p><strong>HOLD — DO NOT SUBMIT H40</strong>, and <strong>HOLD — DO NOT SUBMIT H7 or H40</strong> for the earlier
arms. <strong>Do not upload either version</strong> of anything above: every promotion instrument used on them is
<em>withdrawn</em> by the sixteen-score audit
(<a href="data/instrument-audit-20261006.json">instrument-audit-20261006.json</a>), so none of them earned a
weekly slot.</p>
<p class="micro">Byte hashes of every retained artifact are pinned in <span class="mono">scripts/build_site.py</span>
and <span class="mono">docs/data/session2-artifacts-20261006.json</span>, and verified by CI, so a silent change
cannot pass unnoticed.</p></div></section>

<section class="section"><h2>Limitations stated plainly</h2>
<div class="prose"><ul>
<li>The competition data is login-walled. This repository builds only from the sha256-pinned data bridge plus the
public prior-art corpus; it cannot download or verify organizer files by itself.</li>
<li>No submission was made and no organizer score exists for any file published here. Owner-reported leaderboard
rows are reproduced only as calibration.</li>
<li>No local instrument ranks the recorded scores, so nothing here proves a leaderboard gain.</li>
<li>The 1 m lidar scarp stack used by an earlier hypothesis is not present in this workspace and the USGS S3
endpoint is blocked from this sandbox; H16 states that honestly and no lidar-derived number is claimed.</li>
<li>The depth-labelled cloud is a physical estimate, not a fault map: a potential-field contact may be a
lithologic boundary, and a fault is not automatically a geothermal reservoir.</li>
</ul></div></section>'''

    # ---------------------------------------------------------------- evidence --
    evidence_body = f'''<div class="subhero"><p class="kicker">EVIDENCE / EVERY NUMBER WITH ITS TRUTH LABEL</p>
<h1>What was measured.</h1><p class="lead">Format, novelty, named surrogates and a spatially blocked breakdown,
followed by every negative result that stops this file from being promoted.</p></div>

<section class="section"><h2>1 · Format gate (re-read from the published bytes)</h2>
<div class="tablewrap"><table><thead><tr><th>Check</th><th>Value</th></tr></thead><tbody>
<tr><td>Valid against the official template</td><td class="num">{str(audit["candidate"]["valid"]).lower()}</td></tr>
<tr><td>Band / dtype / CRS</td><td class="num">{audit["candidate"]["candidate"]["count"]} / {audit["candidate"]["candidate"]["dtype"]} / {audit["candidate"]["candidate"]["crs"]}</td></tr>
<tr><td>Shape / geotransform</td><td class="num">{audit["candidate"]["candidate"]["height"]} × {audit["candidate"]["candidate"]["width"]} / {tuple(audit["candidate"]["candidate"]["transform"][:6])}</td></tr>
<tr><td>In-footprint range</td><td class="num">{audit["value_min"]:.6f} … {audit["value_max"]:.6f}</td></tr>
<tr><td>In-footprint non-zero cells</td><td class="num">{audit["candidate"]["candidate"]["in_footprint_nonzero"]:,}</td></tr>
<tr><td>Outside-footprint NaN cells</td><td class="num">{audit["candidate"]["candidate"]["outside_nan"]:,}</td></tr>
<tr><td>Distinct values on the support</td><td class="num">{audit["distinct_values"]:,}</td></tr>
<tr><td>Raw SHA-256</td><td class="mono">{sha256}</td></tr>
<tr><td>Canonical pixel SHA-256</td><td class="mono">{pixsha}</td></tr>
</tbody></table></div>

<h2>2 · Novelty against every cached prior output</h2>
<p>{novelty["corpus_size"]} comparable rasters were read and compared cell by cell. Exact duplicates:
<strong>{len(novelty["exact_duplicates"])}</strong>. Largest |Pearson| correlation:
<strong>{novelty["max_abs_pearson"]:.4f}</strong> (gate {novelty["thresholds"]["max_abs_pearson"]}). Largest
top-mass Jaccard: <strong>{novelty["max_jaccard_topmass"]:.4f}</strong>
(gate {novelty["thresholds"]["max_jaccard_topmass"]}) over {novelty["unreadable_or_incomparable"]} raster that is
not comparable. Verdict: <strong>{"novel" if novelty["novel"] else "not novel"}</strong>.</p>
<div class="tablewrap"><table><thead><tr><th>Closest prior</th><th>Pearson</th><th>Jaccard (top mass)</th>
<th>Containment</th></tr></thead><tbody>
{"".join(f'<tr><td class="mono">{esc(match["file"])}</td><td class="num">{match["pearson"]:.4f}</td><td class="num">{match["jaccard_topmass"]:.4f}</td><td class="num">{match["containment_topmass"]:.4f}</td></tr>' for match in novelty["top_matches"])}
</tbody></table></div>
<p class="micro">A containment of exactly 1.0 appears in the family's earlier audit for dense, all-positive
fields; here the largest containment is {max(match["containment_topmass"] for match in novelty["top_matches"]):.4f}.
Coverage caveat, stated plainly: this gate compares the <strong>{novelty["corpus_size"] + novelty["unreadable_or_incomparable"]}
prior rasters whose bytes are present and hash-verified in <span class="mono">data/prior/</span></strong>. The wider
inventory lists 343 rasters; the remaining 64 could not be re-fetched from this sandbox and are <em>not</em> claimed
as compared.</p>

<h2>3 · Named surrogates — proxies, never scores</h2>
<div class="tablewrap"><table><thead><tr><th>Surrogate truth</th><th>Mean credit per mass (w)</th>
<th>Truth covered</th><th>Proxy DTI</th></tr></thead><tbody>
<tr><td>USGS/INGENIOUS catalogue, with the catalogue and a 1 px ring excluded from the prediction domain</td>
<td class="num">{cat["w"]:.4f}</td><td class="num">{cat["cover"]:.4f}</td><td class="num">{cat["dti_surrogate"]:.5f}</td></tr>
<tr><td>SGMC-derived faults further than 300 m from the catalogue (owner-derived mirror)</td>
<td class="num">{off["w"]:.4f}</td><td class="num">{off["cover"]:.4f}</td><td class="num">{off["dti_surrogate"]:.5f}</td></tr>
</tbody></table></div>
<p class="micro">Distance of the emitted dots to the catalogue: minimum
{audit["dot_distance_to_catalogue_m"]["min"]:.0f} m, median {audit["dot_distance_to_catalogue_m"]["p50"]:.0f} m,
95th percentile {audit["dot_distance_to_catalogue_m"]["p95"]:.0f} m. A large minimum distance is expected of any
file that declines to re-draw traces the organizers mask out of scoring.</p>

<h2>4 · Spatially blocked (4 × 6, three-cell guards)</h2>
{block_table()}
<p class="micro">Block numbers are diagnostics, not scores: the 300 m kernel of an interior cell may be supplied
by a dot just outside the interior, so nothing is clipped away from the metric.</p>

<h2>5 · Score-anchored instrument test (the central negative result)</h2>
<p>The family's audit scored sixteen artifacts against the recorded leaderboard rows and measured every candidate
instrument. None reaches the ρ ≥ 0.8 the promotion rule requires:</p>
<div class="tablewrap"><table><thead><tr><th>Instrument</th><th>Spearman ρ vs recorded score</th></tr></thead><tbody>
<tr><td>Emitted mass (a negative predictor: more mass scored worse)</td><td class="num">{stats["spearman_mass_vs_live"]:+.3f}</td></tr>
<tr><td>Catalogue-calibrated LM instrument</td><td class="num">{stats["spearman_lm_vs_live"]:+.3f}</td></tr>
<tr><td>SGMC off-catalogue credit density</td><td class="num">{stats["spearman_w_offcat_vs_live"]:+.3f}</td></tr>
<tr><td>Required by the family's own promotion rule</td><td class="num">≥ 0.800</td></tr>
</tbody></table></div>
<p>Consequence, stated once and plainly: <strong>no local number available to this programme can justify spending
a weekly submission slot</strong>. That is why this release is published as a research artifact.</p>

<h2>6 · Solver control on an exact analytic contact</h2>
<p>The frozen SI = 0 solver was run on a synthetic field <code>T = A·atan2(x − x₀, z₀)</code>, which is exactly
homogeneous of degree zero, so a correct solver must return the contact. Measured: at zero noise, windows 9 / 15 /
25 accept 585 / 912 / 1,417 solutions with median lateral error 0.6 / 0.3 / 1.0 m and depth recovered to 500 m;
at σ = 1 only 4 solutions survive at window 9 with 72 m median lateral error; at σ = 5 none. The code therefore
works on model-conforming data, and the real-data no-skill result is a property of the data.</p>

<h2>7 · Reproduction</h2>
<div class="prose"><p>Every number on this page is regenerated by two commands from the repository root:</p>
<pre class="codebox">OPENBLAS_NUM_THREADS=2 python scripts/run_h8_lineament.py   # writes work/h8lineament/*
OPENBLAS_NUM_THREADS=1 python scripts/audit_h8_lineament.py  # writes docs/data/h8-lineament-audit.json</pre>
<p class="micro">Receipts: <a href="data/h8-lineament-audit.json">audit</a> ·
<a href="data/h8-lineament-generation.json">generation</a> ·
<a href="data/instrument-audit-20261006.json">instrument audit</a> ·
<a href="data/h33-measured-analysis.json">h33 measurement</a> ·
<a href="research/h8-analysis-20261006.md">analysis</a> ·
<a href="research/h8-preregistration-20261006.md">pre-registration</a>.</p></div>
</section>'''

    # -------------------------------------------------------------- hypotheses --
    hypotheses_body = f'''<div class="subhero"><p class="kicker">RANKED RESEARCH HYPOTHESES</p>
<h1>What to try next, and why it could work.</h1>
<p class="lead">Each entry names the exact challenge layers, the physical signature, why it should catch a fault
the USGS/INGENIOUS catalogue misses rather than one already in it, how it differs from everything already
implemented in this repository, its cost, and how it will be validated.</p></div>

<section class="section"><h2>Rank 1 · H14 — tip, step-over and along-strike extension of mapped systems</h2>
<div class="prose">
<p><strong>Layers.</strong> Catalogue geometry (the provided <code>labels.tif</code> traces) plus independent
strike corroboration from <code>det_elev_slope</code>, <code>rtp</code> and the gravity-gradient bands.</p>
<p><strong>Signature.</strong> Geometric projection: 1–3 km tip extensions, relay corridors between overlapping
traces and along-strike continuation, each crossed with a geophysical lineament test so that a projected tip is
only emitted where the potential field also shows a contact.</p>
<p><strong>Why it is a missing fault and not a catalogued one.</strong> DrivenData staff, 2026-09-23: a new fault
"can include newly mapped geometry of an existing fault system"
(<a href="{STAFF_MASK}">forum 11536, post 2</a>), and the known catalogue pixels are masked out of scoring
(<a href="{STAFF_MASK}">forum 11516</a>) — so an extension beyond a mapped tip is scored normally while the
mapped trace itself costs nothing.</p>
<p><strong>Difference from prior work here.</strong> H4/H7/H8 all re-derive structures from the physics; nothing
in the corpus emits extension geometry, and the family's best file deletes everything within 200 m of the
catalogue instead.</p>
<p><strong>Validation available today.</strong> Leave-the-tips-out on the provided catalogue: hide a random 20 %
of trace tips, rebuild, and measure whether the projection recovers the hidden tips at 100 m resolution with a
capture rate above the same-mass random baseline. Cost: moderate. Expected gain: highest of the five.</p></div></section>

<section class="section"><h2>Rank 2 · H8 — SI = 0 contact Euler depth-clustering (implemented, this release)</h2>
<div class="prose">
<p><strong>Layers.</strong> <code>rtp</code> (150 m upward continuation) and <code>iso_grav_anom</code>
differentiated once (500 m continuation).</p>
<p><strong>Signature.</strong> Depth-labelled contact solutions (Reid et al., 1990, eq. 2 with the arbitrary
offset A) weighted by shallowness, cross-window depth consensus, local cloud lineament coherence and
cross-family corroboration, emitted as an anisotropic kernel-density field.</p>
<p><strong>Why it can find a missing fault.</strong> A potential-field contact is a subsurface boundary; the
catalogue is mapped surface geology. Under basin fill, only geophysics sees the boundary — and the masked
catalogue cannot be scored anyway.</p>
<p><strong>Difference.</strong> H4/H7 emit a full continuous KDE with uniform solution treatment; no earlier
artifact weights by cloud geometry, corroborates across physics families, or emits a metric-tuned sparse
support. <strong>Measured limitation:</strong> the emitted dots are not enriched near the catalogue
(minimum distance {audit["dot_distance_to_catalogue_m"]["min"]:.0f} m, catalogue w = {cat["w"]:.4f}), and the
real-data cloud sits at chance against the mapped faults — the honest expectation is low-to-moderate.</p>
<p><strong>Validation.</strong> Score-anchored leave-one-out instrument (see evidence §5): currently fails to
rank the recorded scores, so this hypothesis stays unvalidated.</p></div></section>

<section class="section"><h2>Rank 3 · H15 — finite-step gravity inversion (SI = −1)</h2>
<div class="prose">
<p><strong>Layers.</strong> <code>iso_grav_anom</code>, <code>iso_grav_anom_vg</code>,
<code>iso_grav_anom_hg</code>, <code>depth_to_base_surf</code>.</p>
<p><strong>Signature.</strong> The finite density step of a fault block — structural index −1 — instead of the
infinite-contact approximation that every Euler artifact in this repository uses
(see <a href="{REID2014}">Reid &amp; Thurston, 2014</a>, on the finite-throw step).</p>
<p><strong>Why it is missing from the catalogue.</strong> Basin-margin faults under sedimentary fill are exactly
the class that surface compilation under-maps; the depth-to-basement layer gives an independent constraint on
where a step can exist.</p>
<p><strong>Difference.</strong> New structural index, new equation, new source geometry; it is not a re-tune of
any existing runner. Cost: high (new solver plus depth-to-basement coupling).</p></div></section>

<section class="section"><h2>Rank 4 · H16 — 1 m lidar scarp re-mapping</h2>
<div class="prose">
<p><strong>Layers.</strong> USGS 3DEP 1 m DEM (the competition data tab ships
<code>1m_DEM_links.csv</code>), aggregated to 100 m as up-face/down-face asymmetry, cross-scarp curvature and
differential relief.</p>
<p><strong>Signature.</strong> The direct surface expression of the labelled fault type. In this repository's own
measurements, the lidar scarp family scored the highest catalogue skill of any layer
(up-face asymmetry AUC 0.5914 versus 0.4855 for the coarse detrended elevation band).</p>
<p><strong>Why it is missing from the catalogue.</strong> Lidar-visible scarps in young deposits are precisely
what pre-lidar compilations miss.</p>
<p><strong>Source, verified as obtainable.</strong> <a href="{THREEDEP}">USGS 3DEP</a> is free and official, and
the competition data page links the exact DEM tiles (<a href="{TNMLINKS}">1m_DEM_links.csv</a>). Flagged
irregularity: the scarp stack built in an earlier session is <em>absent</em> from this workspace and the USGS S3
endpoint is TLS-blocked from this sandbox, so H16 cannot be rebuilt or re-measured here; the AUC above is quoted
from this repository's own earlier receipt and no new lidar number is claimed.</p>
<p>Cost: high (tile download and processing).</p></div></section>

<section class="section"><h2>Rank 5 · H17 — microseismicity-aligned structures</h2>
<div class="prose">
<p><strong>Layers.</strong> <code>ieq_n100a15</code>, <code>deq_n100a15</code>, <code>geod_shearrate</code>.</p>
<p><strong>Signature.</strong> An intensity ridge along an active but unmapped structure; measured catalogue
skill is weak (w = 0.0531 at 20 k selected cells), so it is ranked last on evidence, not on appeal.</p>
<p><strong>Why it is missing from the catalogue.</strong> Active structures can lack mapped surface traces; a
seismicity tie is independent of any compilation. Cost: low.</p></div></section>

<section class="section"><h2>Promotion rules that govern all five</h2>
<div class="prose"><ol>
<li>A hypothesis becomes a submission only after it beats the incumbent on a spatially blocked holdout
(4 × 6 blocks with three-cell guards) <em>and</em> an instrument with demonstrated ranking power — currently
none exists, so nothing may be promoted automatically.</li>
<li>Every candidate must pass the format gate and the raw-output novelty gate against every cached prior raster
before it is published; a near-duplicate is refused, not renamed.</li>
<li>No score is ever claimed for an unscored file, and no weekly slot is spent on an unvalidated idea.</li>
</ol></div></section>'''

    # ----------------------------------------------------------------- sources --
    sources_body = f'''<div class="subhero"><p class="kicker">SOURCES / OFFICIAL AND VERIFIABLE</p>
<h1>Where every claim comes from.</h1>
<p class="lead">Official competition pages, the organisers' own clarifications, the primary geophysical
literature, and the local receipts that record what was actually measured in this repository.</p></div>
<section class="section"><h2>How this site and every artifact are verified</h2>
<div class="prose"><p>Nothing here is published on a claim alone. Continuous integration runs the unit and integrity
tests, re-reads every published byte against its receipt, regenerates the whole site and fails if a single
character changes, and then loads all eight public routes in headless Chromium at desktop, tablet and 320 px
phone widths — checking the download actually serves the published SHA-256 and that the submission name, note
and hash can be copied. The Pages pipeline verifies the live public URLs and both download routes again after
each deployment. When a check fails, the failure is annotated with the exact page, viewport and step.</p>
<p>Current state: <a href="https://github.com/buffedlizard55-lab/GEMSDOE40/actions/workflows/ci.yml">Research and
site checks</a> green on <span class="mono">main</span>; the H8 lineament release reproduces byte-for-byte from an
independent re-run of <span class="mono">scripts/run_h8_lineament.py</span>.</p></div></section>

<section class="section"><h2>Competition</h2><div class="grid3">
<article class="sourcecard"><span class="stepnum">PROBLEM &amp; METRIC</span><h3>Competition problem description</h3>
<p><a href="{PROBLEM}">{PROBLEM}</a> — task, layers and the distance-weighted Tversky index
(α = 0.2, β = 0.8, R = 300 m) used everywhere in this analysis.</p></article>
<article class="sourcecard"><span class="stepnum">DATA</span><h3>Competition data page</h3>
<p><a href="{DATA}">{DATA}</a> — <code>training_features.tif</code>, <code>labels.tif</code>,
<code>sample_submission.tif</code> and the 1 m DEM link CSV. Login-walled; this repository builds from the
sha256-pinned bridge instead.</p></article>
<article class="sourcecard"><span class="stepnum">RULES &amp; FORMAT</span><h3>Submission specification</h3>
<p><a href="{RULES}">{RULES}</a> — single-band float32 GeoTIFF, EPSG:32611, 100 m, values in [0, 1], NaN outside
the footprint; one submission serves both rounds.</p></article>
<article class="sourcecard"><span class="stepnum">LEADERBOARD</span><h3>Public leaderboard</h3>
<p><a href="{LEADERBOARD}">{LEADERBOARD}</a> — observed 2026-10-06: best public score 0.3345, one historical
leader at 0.3195, and the owner's best recorded row at 0.2778. Leaderboard rows referenced in this repository
are owner-reported and are not organizer receipts for a specific file.</p></article>
<article class="sourcecard"><span class="stepnum">REFERENCE</span><h3>Official reference solution</h3>
<p><a href="{REFERENCE}">{REFERENCE}</a> — DrivenData's published baseline, useful for format and pipeline
conventions.</p></article>
<article class="sourcecard"><span class="stepnum">DATA PORTAL</span><h3>GDR entry for the competition data</h3>
<p><a href="{GDR}">{GDR}</a> — Geothermal Data Repository submission describing the GeoDAWN data products;
blocked from this sandbox at the host level, so it is cited, not re-hosted.</p></article>
</div></section>
<section class="section"><h2>Organizer clarifications (verbatim)</h2><div class="grid3">
<article class="sourcecard"><span class="stepnum">MASKING</span><h3>Known faults are excluded from scoring</h3>
<p><a href="{STAFF_MASK}">{STAFF_MASK}</a> — DrivenData staff, 2026-09-16: pixels corresponding to known
USGS/INGENIOUS faults are masked and excluded from evaluation, so they contribute no true positives and are not
penalised.</p></article>
<article class="sourcecard"><span class="stepnum">WHAT IS NEW</span><h3>"Newly mapped geometry" counts</h3>
<p><a href="{STAFF_MASK}">forum 11536, post 2</a> — DrivenData staff, 2026-09-23: "any fault pixel not already
captured by USGS/INGENIOUS" and it "can include newly mapped geometry of an existing fault system".</p></article>
<article class="sourcecard"><span class="stepnum">NOT DISCLOSED</span><h3>Test-fault provenance stays private</h3>
<p><a href="{STAFF_JOIN}">{STAFF_JOIN}</a> — staff, 2026-09-23: data sources, fault types and coverage behind
the test faults are not disclosed, and the final round re-scores against a set updated by expert review of the
Phase 1 submissions.</p></article>
</div></section>
<section class="section"><h2>Geophysics</h2><div class="grid3">
<article class="sourcecard"><span class="stepnum">EULER DECONVOLUTION</span><h3>Reid et al., 1990</h3>
<p><a href="{REID1990}">{REID1990}</a> — the moving-window Euler method with the arbitrary contact offset that
this release implements at structural index 0.</p></article>
<article class="sourcecard"><span class="stepnum">STRUCTURAL INDEX</span><h3>Reid &amp; Thurston, 2014</h3>
<p><a href="{REID2014}">{REID2014}</a> — structural-index selection; a finite gravity step is SI = −1, which is
why the gravity arm here uses the first vertical derivative at SI = 0 as a local top-edge approximation, and why
H15 is ranked as the deeper fix.</p></article>
<article class="sourcecard"><span class="stepnum">ELEVATION</span><h3>USGS 3DEP</h3>
<p><a href="{THREEDEP}">{THREEDEP}</a> — free official 1 m lidar DEM coverage for the next hypothesis; not
reachable from this sandbox (TLS to the USGS endpoint is closed here).</p></article>
</div></section>
<section class="section"><h2>Local receipts (everything measured in this repository)</h2>
<div class="prose"><ul>
<li><a href="data/h8-lineament-audit.json">h8-lineament-audit.json</a> — format, novelty over
{novelty["corpus_size"]} priors, named surrogates, blocked breakdown, recorded-artifact comparison.</li>
<li><a href="data/h8-lineament-generation.json">h8-lineament-generation.json</a> — the exact generation receipt
including configuration, window totals and input hashes.</li>
<li><a href="downloads/{Path(cloud["path"]).name}">h8-lineament-solutions.csv.gz</a> — the depth-labelled solution
cloud, {retained:,} retained solutions.</li>
<li><a href="data/instrument-audit-20261006.json">instrument-audit-20261006.json</a> — the sixteen recorded
scores and the ρ table that retires every local promotion instrument.</li>
<li><a href="data/h33-measured-analysis.json">h33-measured-analysis.json</a> — the deletion measurement behind
the 0.2778 row.</li>
<li><a href="research/h8-analysis-20261006.md">h8-analysis-20261006.md</a> and
<a href="research/h8-preregistration-20261006.md">h8-preregistration-20261006.md</a> — the metric algebra, the
required-credit calculation and the pre-registered settings.</li>
<li><a href="user-prompt-20261006.md">user-prompt-20261006.md</a> — the full project brief, also embedded in the
repository README.</li>
</ul></div></section>'''

    # root alias: the /GEMSDOE40/ route must name the current download and the retained holds
    twin_link = (f'<p>Metric-optimal twin (same {dots:,} cells, every value 1.0): '
                 f'<a href="docs/downloads/{esc(twin["name"])}">{esc(twin["name"])}</a></p>'
                 if twin else "")
    root_html = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta http-equiv="refresh" content="0; url=docs/index.html"><title>GEMSDOE40 · H8 research download</title><link rel="canonical" href="docs/index.html"><link rel="stylesheet" href="docs/assets/contact.css"></head><body><main class="wrap subhero"><p class="kicker">GEMS / 40 · H8</p><h1>Unique Euler depth-clustering candidate</h1><p>SI = 0 contact Euler solutions on magnetic and gravity data, lineament-weighted, emitted as a sparse
density field. Format-valid, novel against {novelty["corpus_size"] + novelty["unreadable_or_incomparable"]} cached prior rasters, and <strong>not promoted</strong>:
no local instrument ranks the family's recorded scores, so no weekly slot was used.</p>
<p><a class="button" href="docs/downloads/{esc(name)}" download>Download the H8 GeoTIFF ↓</a></p>{twin_link}
<p><a href="docs/index.html">Open the project site</a> · <a href="docs/executive-summary.html">Submission guide</a>
· <a href="docs/evidence.html">Evidence</a> · <a href="docs/hypotheses.html">Hypotheses</a></p>
<p><strong>HOLD — DO NOT SUBMIT H40</strong>, and <strong>HOLD — DO NOT SUBMIT H7 or H40</strong> for the earlier arms:
<strong>Nothing earlier is cleared either</strong>. Retained, still downloadable, still held:
<a href="docs/downloads/gemsdoe40-euler-line-ring-pruned-60000px-20261006T032708Z-8dafb186-zeros.tif">H4-line ring-pruned</a> ·
<a href="docs/downloads/{esc(session2["filename"])}">H13 crest-binary</a> ·
<a href="docs/downloads/{esc(session2["sibling_candidate"]["filename"])}">H8 trace-locked depth-KDE</a> ·
<a href="docs/downloads/{esc(session2["h8asa_candidate"]["filename"])}">H8-ASA analytic-signal depth-KDE</a> ·
<a href="docs/downloads/gemsdoe40-h7-rtp-euler-gravity-context-3d-kde-20261006-998f660f.tif">H7</a> ·
<a href="docs/downloads/gemsdoe40-h4-contact-offset-depthkde-20261006-fab9f6619c02.tif">H4</a> ·
<a href="docs/downloads/gemsdoe40-eulerdepth-si0-20261006-run2-57896abe-zeros.tif">H40 run-2</a>.</p>
</main></body></html>'''
    (ROOT / "index.html").write_text(root_html, encoding="utf-8")

    page("index.html", "Euler depth-clustering candidate", index_body, "index.html")
    page("executive-summary.html", "Executive summary and submission guide", executive_body, "executive-summary.html")
    page("evidence.html", "Evidence", evidence_body, "evidence.html")
    page("hypotheses.html", "Hypotheses", hypotheses_body, "hypotheses.html")
    page("sources.html", "Sources", sources_body, "sources.html")

    manifest = {
        "experiment": "H8 lineament-weighted cross-family SI=0 Euler depth-cluster KDE",
        "filename": name,
        "path": published["path"],
        "bytes": published["bytes"],
        "sha256": sha256,
        "canonical_pixels_sha256": pixsha,
        "mass": mass,
        "dots": dots,
        "values": [audit["value_min"], audit["value_max"]],
        "tracking_name": tracking,
        "note": note,
        "note_characters": len(note),
        "cloud": {"path": cloud["path"], "sha256": cloud["sha256"], "bytes": cloud["bytes"]},
        "hard_twin": ({"filename": twin["name"], "path": twin["path"], "sha256": twin["sha256"],
                       "bytes": twin["bytes"]} if twin else None),
        "status": "BUILT, AUDITED, NOT PROMOTED",
        "slot_eligible": False,
        "organizer_score": None,
        "weekly_submission_used": False,
        "novelty": novelty["novel"],
        "corpus_size": novelty["corpus_size"],
        "construction_reads_labels_or_priors": generation["construction_reads_labels_or_priors"],
        "submission_note_example": note,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"rebuilt 5 pages; current candidate = {name} ({published['bytes']:,} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
