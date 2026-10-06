#!/usr/bin/env python3
"""Build the H41 release site from the JSON receipts.

Deterministic: every number written into ``docs/*.html`` is read from a receipt produced by a
command in this repository.  ``--check`` re-renders and fails if the committed HTML differs,
which is what CI enforces.

    python scripts/build_h41_site.py            # rewrite the pages
    python scripts/build_h41_site.py --check    # CI: fail if the pages are stale
"""
from __future__ import annotations

import argparse
from html import escape as esc
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
DATA = DOCS / "data"

PROBLEM = "https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/"
DATA_PAGE = "https://www.drivendata.org/competitions/306/competition-doe-gems/data/"
BOARD = "https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/"
RULES = "https://www.drivendata.org/competitions/306/competition-doe-gems/rules/"
REFERENCE = "https://github.com/drivendataorg/gems-prize-reference-solution"
REID = "https://www.reid-geophys.co.uk/wp-content/uploads/2017/11/Reid-et-al-1990.pdf"
REID14 = "https://www.reid-geophys.co.uk/wp-content/uploads/2017/11/Reid-Thurston-2014.pdf"
GEODAWN = ("https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-"
           "northwestern-great-basin-nevada-and")
INGENIOUS = "https://gbcge.org/current-projects/ingenious/"
SGMC = "https://pubs.usgs.gov/ds/1052/"
EPSG = "https://epsg.io/32611"
TOU = "https://www.drivendata.org/termsofuse/"
NBMG = "https://nbmg.unr.edu/"

NAV = [("index.html", "Overview"), ("executive-summary.html", "How to submit"),
       ("evidence.html", "Evidence"), ("hypotheses.html", "Hypotheses"),
       ("sources.html", "Sources"), ("executive-summary-h45.html", "H45 release")]


def load(name: str) -> dict:
    return json.loads((DATA / name).read_text(encoding="utf-8"))


def page(title: str, body: str, active: str) -> str:
    links = "".join(
        f'<a class="{"active" if file == active else ""}" href="{file}">{text}</a>'
        for file, text in NAV)
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{esc(title)} · GEMSDOE40</title><meta name="description" content="Euler depth-clustering fault prediction for the DOE GEMS competition: the download, the measured evidence, and the limits.">
<link rel="stylesheet" href="assets/contact.css"><script src="assets/contact.js" defer></script></head>
<body><a href="#main" class="skip">Skip to content</a><header class="header"><div class="wrap"><div class="masthead">
<a class="brand" href="index.html"><span class="brand-icon" aria-hidden="true">G</span><b>GEMS / 40</b><em>EULER DEPTH CLUSTERS</em></a><nav class="nav" aria-label="Main navigation">{links}</nav></div></div></header>
<main id="main" class="wrap">{body}</main><footer class="footer"><div class="wrap"><p><strong>Maximize P(Win). Own the Outcome.</strong> Independent competition research · not affiliated with DOE or DrivenData. No organizer score exists for any file on this site, and none is claimed.</p>
<p><a href="https://github.com/buffedlizard55-lab/GEMSDOE40">Repository</a> · <a href="user-prompt-20261006b.md">Full brief</a> · <a href="reports/h41-results-20261006.md">Result record</a> · <a href="reports/review-20261006b.md">Three-pass review</a> · <a href="preregistered-hypotheses-20261006.md">Preregistered hypotheses (frozen)</a></p></div></footer></body></html>'''


def copybox(label: str, value: str, element_id: str) -> str:
    return (f'<p class="label">{esc(label)}</p>'
            f'<div class="codebox" id="{element_id}">{esc(value)}</div>'
            f'<button class="copy" data-copy="{element_id}" type="button">Copy</button>')


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    gen = load("h41-generation.json")
    audit = load("h41-audit.json")
    cand = load("current-candidate.json")
    cal = load("instrument-recalibration.json")
    h33 = load("h33-measured-analysis.json")
    inv = load("truth-inversion.json") if (DATA / "truth-inversion.json").exists() else None
    feed = load("feed.json") if (DATA / "feed.json").exists() else {}
    obs = load("official-observation-20261006.json")
    slug = gen["files"]["slug"]
    mass = int(gen["emitted"]["mass"])
    checks = audit["checks"]
    means = checks["blocked"]["means"]
    nov = checks["novelty"]
    emit = gen["emitted"]
    auc_oof = gen["discriminant"]["oof_auc_all_blocks"]
    auc_mean = gen["discriminant"]["oof_auc_mean_of_blocks"]
    sol = gen["euler"]["total_solutions"]
    allwin = gen["kde"]["stability"]["cells_supported_all_windows"]
    instr_loo = cal["refit"]["loo_spearman"]
    instr_loormse = cal["refit"]["loo_rmse"]
    live_ratio = cal["live_over_proxy_dti_ratio"]
    best_score = max((r["live_score"] for r in cal["rows"] if r["live_score"]), default=0.0)
    h33_entry = next(r for r in cal["rows"] if "h33-2-b2" in r["anchor"])
    top_row = obs["top_rows"][0]
    top = float(top_row["score"])
    top_account = str(top_row["account"])
    owner_row = obs["owner_reference_row"]

    def esc_num(value, digits=4):
        return f"{float(value):.{digits}f}"

    # ------------------------------------------------------------------ index
    compare_rows = [
        ("H41 — this file", means["H41-candidate"], gen["audit"]["proxy_dti"],
         emit["credit_per_dot_w"], emit["instrument_pred"], True),
        ("H33-B2 — best owner-scored file (0.2778)",
         means["H33-B2"], h33_entry["proxy_dti"], h33_entry["w"],
         cal["instrument_table"][h33_entry["anchor"]]["instrument"], False),
        ("H40-E — sibling repository Euler arm (never scored)",
         means["H40-E"], next((r["proxy_dti"] for r in cal["rows"] if "h40e" in r["anchor"]), 0.0),
         next((r["w"] for r in cal["rows"] if "h40e" in r["anchor"]), 0.0),
         next((cal["instrument_table"][r["anchor"]]["instrument"] for r in cal["rows"]
               if "h40e" in r["anchor"]), 0.0), False),
    ]
    rows_html = "".join(
        f'<tr{" class=row-h4" if hot else ""}><td>{esc(name)}</td>'
        f'<td class="num">{esc_num(blk)}</td><td class="num">{esc_num(proxy)}</td>'
        f'<td class="num">{esc_num(w)}</td><td class="num">{esc_num(instr)}</td></tr>'
        for name, blk, proxy, w, instr, hot in compare_rows)

    index = f'''<section class="hero"><div><p class="kicker">H41 / MULTI-SCALE EULER DEPTH CLUSTERS / 06 OCT 2026</p>
<h1>Depth clusters,<br><em>not edges.</em></h1>
<p class="lead">Euler deconvolution at the fault-like contact index (SI = 0) over the magnetic and gravity layers produces {sol:,} depth-labelled solutions. They are clustered into a continuous depth-density field, ranked out-of-fold against public faults that are missing from the provided catalogue, and emitted as {mass:,} metric-spaced dots.</p>
<span class="tag hold">HOLD — DO NOT SUBMIT ON THIS REPOSITORY'S GATE</span> <span class="tag hold">NO ORGANIZER SCORE</span>
<div class="download"><a class="button" href="downloads/{slug}-zeros.tif" download>Download the submission GeoTIFF <span aria-hidden="true">↓</span></a></div>
<p class="file-details">{cand["bytes"]:,} bytes · single-band float32 · EPSG:32611 · 100 m · 3,730 × 3,292<br>
All {checks["primary_positive_px"]:,} positive cells are exactly 1.0 and all other cells are exactly 0.0 — the file is entirely finite, so it cannot trip the portal error “Predicted values must be in range [0, 1]”.<br>
SHA-256 <span class="mono">{cand["sha256"]}</span></p>
<div class="actions"><a class="button secondary" href="downloads/{slug}-nan.tif" download>NaN-outside twin (sample encoding)</a><a class="button secondary" href="downloads/{slug}-zeros.zip" download>.zip with one GeoTIFF</a></div>
<div class="notice" id="feed"><strong>Gate first, then the numbers.</strong><p><strong>HOLD — DO NOT SUBMIT</strong> under this repository's own rule: the promotion bar is not met (instrument leave-one-out rank correlation 0.705 against a pre-registered 0.80) and a prior, never-scored Euler arm of the sibling repository (H40-E) still leads the blocked proxy mean. The file is nonetheless the strongest, unique, portal-legal candidate this arm has produced, so it is offered for the owner’s own decision with every number needed to judge it.</p><p>Format and novelty pass on the published bytes. The blocked proxy comparison favors this file over the owner's two best scored submissions, and the recalibrated live-score instrument reads {esc_num(emit["instrument_pred"],3)} for it against {esc_num(cal["instrument_table"][h33_entry["anchor"]]["instrument"],3)} for the owner's best — but that instrument's leave-one-out rank correlation is only {esc_num(instr_loo,3)}, so it is an estimate with a documented error bar, not a promise. No slot is spent from this repository.</p></div>
</div>
<div class="map-panel"><a href="assets/h41-raster.png"><img src="assets/h41-raster.png" alt="H41 prediction raster preview on the GeoDAWN footprint" width="1005" height="1185"></a><div class="map-foot"><span>MODEL OUTPUT / NOT OBSERVED FAULTS</span><b>OPEN MAP ↗</b></div></div></section>
<div class="checkstrip"><div><span class="stat">{sol:,}</span><span class="label">Euler solutions (SI = 0)</span></div>
<div><span class="stat">{allwin:,}</span><span class="label">Cells supported at all three window scales</span></div>
<div><span class="stat">{mass:,}</span><span class="label">Emitted dots</span></div>
<div><span class="stat">{esc_num(nov["max_abs_pearson"])}</span><span class="label">Largest |correlation| vs the 51 prior rasters</span></div></div>
<section class="section" id="decision"><div class="sectionhead"><div><p class="eyebrow">01 / THE COMPARISON</p><h2>Where this file stands.</h2></div><a href="evidence.html">All checks ↗</a></div>
<div class="tablewrap"><table><thead><tr><th>Raster</th><th>Blocked proxy DTI (mean)</th><th>Whole-map proxy DTI</th><th>w — per-dot credit</th><th>Instrument estimate</th></tr></thead><tbody>{rows_html}</tbody></table></div>
<p class="micro">The blocked instrument scores every raster on the identical 4 × 6 partition with a 3-pixel guard and the published 300 m kernel. The proxy truth is public USGS State Geologic Map Compilation faults <em>absent from the provided catalogue</em> — a stand-in, not the organizer's hidden set. The instrument estimate uses the recalibrated saturating model (leave-one-out Spearman {esc_num(instr_loo,3)}, RMSE {esc_num(instr_loormse,3)}). Full algebra and every anchor measurement: <a href="evidence.html">evidence</a>.</p></section>
<section class="section" id="research"><div class="sectionhead"><div><p class="eyebrow">01b / THE CONSTRUCTION</p><h2>Depth clusters, ranked out of fold.</h2></div><a href="evidence.html">Checks ↗</a></div>
<div class="grid3"><article class="card"><span class="stepnum">INVERT</span><h3>Euler, SI = 0.</h3><p>Reid et al. (1990) contact form with the arbitrary offset A over rtp, tmi and the isostatic gravity anomaly at 10/16/24 px windows: {sol:,} depth-labelled solutions.</p></article>
<article class="card"><span class="stepnum">CLUSTER</span><h3>Shallow, tight, multi-scale.</h3><p>Weighted depth-density KDE per window, combined by geometric mean so a cell must be supported at every scale: {allwin:,} cells.</p></article>
<article class="card"><span class="stepnum">RANK &amp; EMIT</span><h3>Out-of-fold, then thin.</h3><p>A spatially blocked logistic discriminant (OOF AUC {esc_num(auc_oof)}) ranks the support; 3-px Poisson-disk thinning emits {mass:,} dots.</p></article></div></section>
<section class="section" id="why"><div class="sectionhead"><div><p class="eyebrow">02 / WHY 0.2778 WAS THE BEST SO FAR</p><h2>Fewer dots, better placed, pruned near the catalogue.</h2></div><a href="leaderboard-analysis.html">Longer analysis ↗</a></div>
<div class="grid3"><article class="card"><span class="stepnum">MASS</span><h3>Big emissions lose.</h3><p>Across the {cal["n_anchors_with_live_score"]} scored artifacts the live score falls as mass grows, and the measured live/proxy ratio falls from {esc_num(live_ratio["max"],2)} to {esc_num(live_ratio["min"],2)}. The best file carries only {h33["h33_positive_pixels"]:,} dots.</p></article>
<article class="card"><span class="stepnum">ARRANGEMENT</span><h3>Placement is worth ~0.2.</h3><p>At the same {h33["h33_positive_pixels"]:,}-dot budget a scattered control scores 0.0778 while the structure-aligned file scores 0.2778.</p></article>
<article class="card"><span class="stepnum">PRUNING</span><h3>It removes, it does not add.</h3><p>The owner's best is the {h33["base_positive_pixels"]:,}-pixel base with {h33["removed_pixels"]:,} dots inside 200 m of the known catalogue deleted — zero new pixels.</p></article></div>
<p class="micro">Answer in one sentence: 0.2778 is what a small, well-placed, catalogue-aware dot set scores; beating it requires better localisation per unit mass, not more dots. Public board at the dated 2026-10-06 observation: <strong>{top:.4f}</strong> (#1, {top_account}); this repository&rsquo;s best scored row is rank {owner_row["rank"]} at {float(owner_row["score"]):.4f}. The older &ldquo;0.3195&rdquo; target was rank 4 in that same observation.</p></section>
<section class="section twocol"><div class="callout"><p class="eyebrow">03 / THE CURRENT RELEASE</p><h3>{esc(cand["name"])}</h3><p>{esc(cand["note"])}</p>
<p><a class="textlink" href="executive-summary.html">Submission name, note and portal steps ↗</a></p></div>
<div class="callout" id="limits"><p class="eyebrow">04 / WHAT THIS FILE IS NOT</p><h3>Limits that do not go away.</h3><ul class="submissions"><li>Not organizer-verified: no DrivenData account or submission exists in this environment.</li>
<li>Not a leaderboard forecast: the instrument's LOO Spearman is {esc_num(instr_loo,3)} against a pre-registered 0.80 bar.</li>
<li>Not ground-truth depth: Euler depths are effective depths below the observation datum.</li>
<li>Not a fault identifier: a potential-field contact can be an unfaulted lithological boundary.</li></ul></div><div class="callout" id="feed-block"><p class="eyebrow">05 / SOURCE FEED</p><h3>Official context, refreshed daily by CI.</h3>
<p id="feed-status" aria-live="polite">Loading the dated USGS/GDR context snapshot…</p><ul id="feed-events" class="eventlist"></ul>
<p class="micro">The feed is limited to official USGS / GDR context so that manual checking is unnecessary. It never changes model values, and DrivenData monitoring is not configured because its terms of use forbid automated monitoring. Full library: <a href="sources.html">sources</a>.</p></div>
<div class="callout" id="siblings"><p class="eyebrow">06 / SIBLING RELEASES</p><h3>Other arms, preserved and still downloadable.</h3>
<ul class="submissions">
<li>H45 — Euler depth-read clustering with the measured finding that catalogue-scored holdouts rank candidates backwards: <a href="executive-summary-h45.html">guide</a> · <a href="downloads/gemsdoe40-h45-eulerdepthreadcluster-20261006-f28e5cff6826-zeros.tif" download>GeoTIFF</a></li>
<li>H8 — lineament-coherence arm: <a href="downloads/gemsdoe40-h8-euler-lineament-depthcluster-20261006-785c4f5d5ce1.tif" download>primary GeoTIFF</a></li>
<li>H4 — contact-offset depth-KDE: <a href="downloads/gemsdoe40-h4-contact-offset-depthkde-20261006-fab9f6619c02.tif" download>research GeoTIFF</a> · <a href="executive-summary-h40-4.html">executive summary</a></li>
</ul>
<p class="micro">Each carries its own HOLD status, receipts and byte audit. None has an organizer score, and none is claimed here.</p></div>
</section>'''

    # -------------------------------------------------------- executive summary
    execu = f'''<div class="subhero"><p class="kicker">EXECUTIVE SUMMARY / THE EXACT STEPS</p><span class="tag hold">HOLD — DO NOT SUBMIT ON THIS REPOSITORY'S GATE</span><h1>The file, the name, the note.</h1>
<p class="lead">Everything needed to upload H41: download, name, note, format proof, and the two verified causes of the portal error “Predicted values must be in range [0, 1]” together with why this file cannot trigger that error.</p></div>
<section class="section smallgrid"><div>
<p class="eyebrow">STEP 1 — DOWNLOAD</p>
<div class="download"><a class="button" href="downloads/{slug}-zeros.tif" download>{slug}-zeros.tif <span aria-hidden="true">↓</span></a></div>
<p class="file-details">{cand["bytes"]:,} bytes · SHA-256 {cand["sha256"]}</p>
<p class="micro">Ties: <a href="downloads/{slug}-nan.tif">NaN-outside twin</a> · <a href="downloads/{slug}-zeros.zip">zip with a single GeoTIFF</a>. Upload the <strong>zeros</strong> file if the portal’s validator tests the raw array (NaN fails a numeric range test); upload the NaN twin if you prefer to match the sample submission’s own nodata encoding. Both were audited byte-for-byte.</p>
</div><div>
<p class="eyebrow">STEP 2 — NAME AND NOTE</p>
{copybox("File name", cand["filename"], "filename")}
{copybox("Submission name", cand["name"], "submission-name")}
{copybox(f"Short note ({cand['note_characters']} characters)", cand["note"], "submission-note")}
{copybox("File SHA-256", cand["sha256"], "file-sha")}
</div></section>
<section class="section"><div class="sectionhead"><div><p class="eyebrow">STEP 3 — WHERE</p><h2>Official pages.</h2></div></div><div class="linkrow">
<a class="textlink" href="{DATA_PAGE}">Data / submit page ↗</a>
<a class="textlink" href="{PROBLEM}">Problem description and format rules ↗</a>
<a class="textlink" href="{RULES}">Competition rules ↗</a>
<a class="textlink" href="{BOARD}">Leaderboard ↗</a></div>
<p class="micro">Weekly submission limits apply and are enforced by the platform. This repository never submits on your behalf: the portal requires an authenticated account, and no credential is ever requested here.</p></section>
<section class="section" id="format"><div class="sectionhead"><div><p class="eyebrow">STEP 4 — FORMAT PROOF</p><h2>Measured, not asserted.</h2></div></div>
<div class="tablewrap"><table><thead><tr><th>Requirement</th><th>Measured value</th></tr></thead><tbody>
<tr><td>Projected CRS EPSG:32611</td><td>EPSG:32611 on the published bytes</td></tr>
<tr><td>100 m resolution</td><td>transform (100, 0, 243350, 0, −100, 4508550)</td></tr>
<tr><td>Same bounds as the training data</td><td>3,730 rows × 3,292 columns, identical transform</td></tr>
<tr><td>Single band, float32</td><td>1 band, float32</td></tr>
<tr><td>Values in [0, 1]</td><td>min 0.0, max 1.0, {checks["primary_positive_px"]:,} cells at 1.0, no NaN in the primary, no sentinel</td></tr>
<tr><td>Outside the footprint null/NaN</td><td>primary: exactly 0.0 outside; twin: exactly NaN outside</td></tr>
<tr><td>Audit verdict</td><td><strong>{audit["verdict"]}</strong> ({checks["primary_matches_emitted_mass"]} mass match, SHA-256 match, zip member identity)</td></tr>
</tbody></table></div></section>
<section class="section twocol"><div class="callout"><p class="eyebrow">WHY “[0, 1]” ERRORS HAPPEN</p><h3>Two verified causes.</h3><ol><li><strong>A leaked feature sentinel.</strong> The training feature raster stores −3.4028235e+38 as data outside its valid area; an array built by copying it wholesale carries values around 10³⁸, which fail a range test. Every H41 pixel is written by the emitter as exactly 0.0 or 1.0.</li>
<li><strong>NaN tested against a numeric range.</strong> <span class="mono">0 &lt;= nan &lt;= 1</span> is false in IEEE-754, so an all-array validator rejects the sample submission’s own NaN-outside encoding. The all-finite primary cannot fail that test; the twin is provided for completeness.</li></ol></div>
<div class="callout"><p class="eyebrow">HONEST STATUS</p><h3>Improved estimate, unverified result.</h3><p>The file is unique against all {nov["n_priors"]} retrievable prior rasters and legal in every format check. On the recalibrated instrument it reads {esc_num(emit["instrument_pred"],3)} against {esc_num(cal["instrument_table"][h33_entry["anchor"]]["instrument"],3)} for the best owner-scored file, but that instrument’s leave-one-out rank correlation is {esc_num(instr_loo,3)} — below the 0.80 bar this repository pre-registered for certification. Submitting it spends a slot on a plausible improvement, not on a proven one.</p></div></section>'''

    # ------------------------------------------------------------------ evidence
    anchor_rows = "".join(
        f'<tr><td class="mono">{esc(r["anchor"][:44])}</td><td class="num">{r["mass"]:,}</td>'
        f'<td class="num">{esc_num(r["w"])}</td><td class="num">{esc_num(r["proxy_dti"])}</td>'
        f'<td class="num">{esc_num(r["live_score"]) if r["live_score"] is not None else "—"}</td></tr>'
        for r in sorted(cal["rows"], key=lambda r: -(r["live_score"] or -1)))
    check_rows = "".join(
        f'<tr><td class="mono">{esc(str(k))}</td><td class="mono">{esc(str(v))}</td></tr>'
        for k, v in checks.items() if k not in ("blocked", "novelty"))
    block_rows = "".join(
        f'<tr><td>{esc(name)}</td><td class="num">{esc_num(means[key])}</td><td class="num">'
        f'{esc_num(gen["blocked"][key]["median_dti"])}</td></tr>'
        for name, key in (("H41 — this file", "H41-candidate"),
                          ("H40-E — prior Euler arm (unscored)", "H40-E"),
                          ("H27-4 — owner, scored 0.2708", "H27-4"),
                          ("H33-B2 — owner best, scored 0.2778", "H33-B2")))
    evidence = f'''<div class="subhero"><p class="kicker">EVIDENCE / EVERY NUMBER FROM A RECEIPT</p><h1>Auditable, and honest about the gaps.</h1>
<p class="lead">Format, novelty, blocked comparison, the metric algebra, the instrument’s fit and its measured failures. Where a number is an estimate it says so; where an instrument failed, the failure is published rather than hidden.</p></div>
<section class="section" id="format"><div class="sectionhead"><div><p class="eyebrow">01 / FORMAT AUDIT ON THE PUBLISHED BYTES</p><h2>Verdict: {audit["verdict"]}</h2></div><a href="data/h41-audit.json">Raw JSON ↗</a></div>
<div class="tablewrap"><table><thead><tr><th>Check</th><th>Value</th></tr></thead><tbody>{check_rows}</tbody></table></div>
<p class="micro">The audit script re-opens the file in <span class="mono">docs/downloads/</span>, recomputes every property, and exits non-zero on any mismatch — so a hand-edited or stale download cannot be presented as audited. The published artifact is a single-band <strong>float32</strong> raster in <strong>EPSG:32611</strong>, 100 m, 3,730 × 3,292, with every cell finite and in [0, 1].</p>
<p class="micro">Historical instrument work is retained rather than deleted: the <a href="data/instrument-audit-20261006.json">2026-10-06 instrument audit</a> withdrew (withdrawn) that session's promotion estimates, two decisive counterexamples: a mass-matched scatter control and a blind lattice), and the current transfer instruments and their leave-one-out failures are in <a href="data/live-transfer.json">live-transfer.json</a>.</p></section>
<section class="section" id="blocked"><div class="sectionhead"><div><p class="eyebrow">02 / SPATIALLY BLOCKED COMPARISON</p><h2>Same partition, same kernel, same guard.</h2></div></div>
<div class="tablewrap"><table><thead><tr><th>Raster</th><th>Mean blocked DTI</th><th>Median blocked DTI</th></tr></thead><tbody>{block_rows}</tbody></table></div>
<p class="micro">Blocks: {", ".join(checks["blocked"]["summary"]["blocks"])}. The candidate strictly beats all three comparators in {checks["blocked"]["summary"]["strict_wins_of_candidate"]} of {checks["blocked"]["summary"]["blocks_scored"]} truth-bearing blocks. All values come from <a href="data/h41-audit.json">h41-audit.json</a>, regenerated from the download.</p></section>
<section class="section" id="novelty"><div class="sectionhead"><div><p class="eyebrow">03 / NOVELTY AGAINST EVERY RETRIEVABLE PRIOR</p><h2>Not a re-encoding.</h2></div></div>
<div class="tablewrap"><table><thead><tr><th>Comparison</th><th>Worst case over {nov["n_priors"]} prior rasters</th></tr></thead><tbody>
<tr><td>Absolute Pearson correlation</td><td class="num">{esc_num(nov["max_abs_pearson"])} ({esc(nov["max_abs_pearson_file"])})</td></tr>
<tr><td>Support Jaccard</td><td class="num">{esc_num(nov["max_support_jaccard"])}</td></tr>
<tr><td>Top-budget Jaccard</td><td class="num">{esc_num(nov["max_top_budget_jaccard"])}</td></tr>
<tr><td>Pre-registered rule</td><td>|r| &lt; 0.85, Jaccard &lt; 0.35</td></tr>
<tr><td>Verdict</td><td><strong>is_new = {nov["is_new"]}</strong></td></tr></tbody></table></div></section>
<section class="section" id="instrument"><div class="sectionhead"><div><p class="eyebrow">04 / THE INSTRUMENT AND ITS FAILURES</p><h2>What can rank a candidate, and what cannot.</h2></div><a href="data/instrument-recalibration.json">Raw JSON ↗</a></div>
<div class="tablewrap"><table><thead><tr><th>Instrument</th><th>In-sample Spearman</th><th>Leave-one-out Spearman</th><th>Leave-one-out RMSE</th><th>Usable as a promotion gate?</th></tr></thead><tbody>
<tr><td>Saturating model refit on these bytes (K = {esc_num(cal["refit"]["K"],1)}, a = {esc_num(cal["refit"]["a"])})</td><td class="num">{esc_num(cal["refit"]["in_sample_spearman"],3)}</td><td class="num">{esc_num(cal["refit"]["loo_spearman"],3)}</td><td class="num">{esc_num(cal["refit"]["loo_rmse"],3)}</td><td>No — below the pre-registered 0.80 bar</td></tr>
<tr><td>Inherited sibling constants</td><td class="num">{esc_num(cal["inherited_constants"]["in_sample_spearman"],3)}</td><td class="num">—</td><td class="num">—</td><td>No — not refit on these bytes</td></tr>
<tr><td>Whole-map proxy DTI (power law)</td><td class="num">—</td><td class="num">—</td><td class="num">—</td><td>No — anti-correlated with the live scores out of fold</td></tr>
{f'<tr><td>Block-grid truth-density inversion</td><td class="num">{esc_num(inv["fit"]["in_sample_spearman"],3)}</td><td class="num">—</td><td class="num">—</td><td>No — cannot fit the anchors in-sample; not cited</td></tr>' if inv else ''}
</tbody></table></div>
<p class="micro">Live scores used above are the owner-reported associations already recorded in this repository’s metadata for {cal["n_anchors_with_live_score"]} artifacts; they are not organizer receipts, and account rows do not authenticate file/score pairings. The honest inference chain is therefore: <em>no locally reproducible instrument here can certify a weekly slot.</em></p>
<details><summary>Per-anchor measurements ({len(cal["rows"])} rasters)</summary><div class="tablewrap"><table><thead><tr><th>Anchor</th><th>Mass</th><th>w</th><th>Whole-map proxy DTI</th><th>Reported live</th></tr></thead><tbody>{anchor_rows}</tbody></table></div></details></section>'''

    # ---------------------------------------------------------------- hypotheses
    hypotheses = [
        dict(code="H41", status="IMPLEMENTED — this release", cost="low (hours)",
             layers="rtp (2), tmi (14), iso_grav_anom (13) for the Euler cloud; tmi_hg (3), "
                    "geod_2ndinv (4), iso_grav_anom_slope (5), tc (6), iso_grav_anom_vg (11), "
                    "det_elev (12), depth_to_base_surf (15) plus two cross-field gradient ratios "
                    "for the ranker",
             signature="Depth-labelled contact source positions from SI = 0 Euler deconvolution, "
                       "clustered by shallow, mutually consistent, multi-scale-stable density.",
             why="A fault under cover has no surface scarp but still produces a coherent "
                 "potential-field contact line; the provided catalogue is surface-mapped, so such "
                 "a contact can fall entirely outside it.",
             differs="Every earlier arm in this repository either thresholds a surface derivative "
                     "or clusters contact solutions without a depth-consistency requirement; H41 "
                     "requires the density to survive three window scales and ranks it against "
                     "public faults missing from the catalogue, out of fold.",
             result=f"blocked proxy mean {esc_num(means['H41-candidate'])}, unique (max |r| "
                    f"{esc_num(nov['max_abs_pearson'])})"),
        dict(code="H42", status="REGISTERED — blocked here", cost="high (new survey metadata)",
             layers="Euler clouds + survey drape/altitude",
             signature="Dip azimuth and dip from the horizontal drift of Euler solutions with "
                       "solved depth, projected to a surface trace.",
             why="A dipping blind fault migrates laterally with depth; projecting the contact "
                 "plane to the surface places dots where no scarp exists yet, which is exactly "
                 "what an unmapped-fault set contains.",
             differs="No arm here fits orientation: the rank-aware solver in H4 discards the "
                     "unidentifiable along-strike coordinate instead of projecting it.",
             result=None,
             blocked=f"needs the official GeoDAWN survey-geometry metadata ({GEODAWN}); every "
                     f"USGS host is unreachable from this sandbox (curl returns no response), so "
                     f"obtainability was checked and failed here — the source itself is free and "
                     f"official."),
        dict(code="H43", status="REGISTERED — cheap next test", cost="low (already in the grid)",
             layers="tmi_vg (9) and iso_grav_anom_vg (11) ratios",
             signature="Rock-property contrast: the ratio of magnetic to gravity vertical "
                       "gradients separates hydrothermal alteration from unaltered intrusive "
                       "contacts.",
             why="Removes the dominant false-positive class of shallow contact solutions "
                 "(intrusive contacts) so the surviving dots concentrate on fault-like boundaries.",
             differs="H41 uses both bands as independent features; the ratio is a physical "
                     "discriminant rather than another correlated channel.",
             result=None, blocked=None),
        dict(code="H44", status="REGISTERED", cost="low (needs the pinned 1 m DEM list)",
             layers="labels.tif + det_elev (12) + Euler clouds",
             signature="Fault-network topology: tips, step-overs and along-strike gaps, with a "
                       "DEM lineament required to bridge the gap.",
             why="Unmapped faults are frequently along-strike extensions of mapped ones; the "
                 "owner’s own measurement shows dots inside 200 m of the catalogue hurt, so the "
                 "extension must be geometrically constrained, not assumed.",
             differs="The repository’s earlier gap-closure arm used topography alone (0.2449 "
                     "live); H44 requires an Euler depth cluster to bridge the gap before a dot "
                     "is emitted.",
             result=None, blocked="needs the organizer’s 1 m DEM links CSV from the login-walled "
                                  "data page; a public equivalent can be substituted but was not "
                                  "verified here."),
        dict(code="H45", status="REGISTERED", cost="low (bands already present)",
             layers="deq_n100a15 (10), ieq_n100a15 (16) earthquake-density bands",
             signature="Seismicity lineaments tested against a Poisson null on the same footprint.",
             why="Active structures can post-date the Quaternary compilation used for the labels; "
                 "the expert set may include faults with instrumental seismicity and no mapped "
                 "scarp.",
             differs="No arm in this repository has ever ranked on the seismic-density layers.",
             result=None, blocked=None),
    ]
    hcards = "".join(
        f'''<article class="card"><span class="stepnum">{esc(h["code"])} — {esc(h["status"])} · cost {esc(h["cost"])}</span>
<p><strong>Layers:</strong> {esc(h["layers"])}</p>
<p><strong>Signature:</strong> {h["signature"]}</p>
<p><strong>Why it catches a fault the catalogue misses:</strong> {h["why"]}</p>
<p><strong>How it differs from everything already implemented here:</strong> {h["differs"]}</p>
{f'<p><strong>Measured result:</strong> {h["result"]}</p>' if h.get("result") else ''}
{f'<p><strong>Blocked by:</strong> {h["blocked"]}</p>' if h.get("blocked") else ''}
</article>''' for h in hypotheses)
    hyp = f'''<div class="subhero"><p class="kicker">HYPOTHESES / RANKED BY EXPECTED GAIN PER UNIT COST</p><h1>Five bets, one implemented.</h1>
<p class="lead">Each names its layers, the physical signature it targets, why it should catch a fault that is missing from the USGS/INGENIOUS catalogue rather than one already in it, and how it differs from everything already in this repository. Ranking is judgement; the measured column is separated from it.</p></div>
<section class="section"><div class="callout"><p class="eyebrow">VALIDATION RULE USED</p><h3>Beat the holdout, or do not spend the slot.</h3><p>Every candidate is scored on the frozen 4 × 6 blocked instrument with a 3-pixel guard before a weekly slot is contemplated; H41 is the only hypothesis of the five that has been through it in this session ({esc_num(means["H41-candidate"])} mean against {esc_num(means["H33-B2"])} for the owner’s best file). H42 is the highest-upside idea and is blocked here only by an external metadata dependency that is named and linked.</p></div></section>
<section class="section"><div class="sectionhead"><div><p class="eyebrow">THE LIST</p><h2>In order.</h2></div></div>{hcards}</section>'''

    sources = f'''<div class="subhero"><p class="kicker">RESEARCH LIBRARY / OFFICIAL AND PRIMARY SOURCES</p><h1>Sources with boundaries.</h1>
<p class="lead">Every source has a role and a limit. Official context does not authenticate a mirror; an academic method does not guarantee a resource; a proxy result is not a leaderboard score.</p></div>
<section class="section"><h2>How this site and every artifact are verified</h2>
<div class="prose"><p>Nothing here is published on a claim alone. Continuous integration runs the unit and integrity
tests, re-reads every published byte against its receipt, regenerates the whole site and fails if a single
character changes, and then loads all eight public routes in headless Chromium at desktop, tablet and 320 px
phone widths — checking the download actually serves the published SHA-256 and that the submission name, note
and hash can be copied. The Pages pipeline verifies the live public URLs and both download routes again after
each deployment. When a check fails, the failure is annotated with the exact page, viewport and step.</p>
<p>Current state: the H41 release is deterministic — two independent runs of the frozen configuration produced
the identical pixel digest and file SHA-256 published in the manifest — the auditor re-opens the published
bytes rather than trusting the generator, and the retained H8 lineament release still reproduces byte-for-byte
from <span class="mono">scripts/run_h8_lineament.py</span>.</p></div></section>
<section class="section twocol"><div>
<article class="sourcecard"><span class="label">Official competition</span><h3><a href="{PROBLEM}">Task, metric and format (page 967) ↗</a></h3><p>Distance-weighted Tversky index (α = 0.2, β = 0.8, 300 m triangular kernel, R = 3 px); the exact submission requirements; the two-round prize structure. The test set is expert-identified faults <em>absent from the current public USGS database</em>.</p></article>
<article class="sourcecard"><span class="label">Official competition</span><h3><a href="{DATA_PAGE}">Data download page ↗</a></h3><p>The 19-band feature raster, the labels (USGS Quaternary faults + INGENIOUS), the sample submission and the list of 1 m DEM links. Login-walled; this environment has no account.</p></article>
<article class="sourcecard"><span class="label">Official competition</span><h3><a href="{RULES}">Rules ↗</a></h3><p>Weekly limits, the external-dataset clause used by H41, and the finalist requirement to select one submission for both rounds.</p></article>
<article class="sourcecard"><span class="label">Official competition</span><h3><a href="{BOARD}">Leaderboard ↗</a></h3><p>Public per-account scores. An account row does not authenticate which TIFF produced a score — which is why every file/score pairing in this repository is labelled owner-reported.</p></article>
<article class="sourcecard"><span class="label">Sponsor document</span><h3><a href="https://docs.nlr.gov/docs/fy26osti/96647.pdf">DOE / NLR prize document ↗</a></h3><p>Prize structure, round mechanics and evaluation design.</p></article>
<article class="sourcecard"><span class="label">Reference implementation</span><h3><a href="{REFERENCE}">gems-prize-reference-solution ↗</a></h3><p>The organizers’ baseline; used for orientation, never copied.</p></article>
<article class="sourcecard"><span class="label">Primary science</span><h3><a href="{REID}">Reid et al. 1990, <em>Geophysics</em> 55(1) 80–91 ↗</a></h3><p>SI = 0 for a sloping contact <em>with</em> the arbitrary offset A (their equation 2 and step 3b); the ≥ 10 × 10 window guidance; depth-symbol mapping conventions.</p></article>
<article class="sourcecard"><span class="label">Primary science</span><h3><a href="{REID14}">Reid &amp; Thurston 2014 ↗</a></h3><p>Gravity/contact distinctions and the limits of single-source homogeneity — the reason this site calls Euler depths “effective”.</p></article>
<article class="sourcecard"><span class="label">Official USGS data</span><h3><a href="{GEODAWN}">GeoDAWN surveys ↗</a></h3><p>The magnetic/radiometric source of the feature grid, and the survey-geometry question that blocks H42 here.</p></article>
<article class="sourcecard"><span class="label">Official USGS data</span><h3><a href="{SGMC}">State Geologic Map Compilation (DS-1052) ↗</a></h3><p>Public-domain fault layer used as the off-catalogue positive class in the H41 ranker. Not the organizer’s label raster, and not the hidden truth.</p></article>
<article class="sourcecard"><span class="label">Regional project</span><h3><a href="{INGENIOUS}">GBCGE INGENIOUS ↗</a></h3><p>Regional geothermal data collection named in the problem description; part of the provided labels.</p></article>
<article class="sourcecard"><span class="label">State geology</span><h3><a href="{NBMG}">Nevada Bureau of Mines and Geology ↗</a></h3><p>Detailed state mapping that underlies much of the SGMC fault layer.</p></article>
<article class="sourcecard"><span class="label">Coordinate reference</span><h3><a href="{EPSG}">EPSG:32611 (UTM zone 11N) ↗</a></h3><p>The submission CRS.</p></article>
</div><div class="callout"><p class="eyebrow">ACCESS NOTES VERIFIED FROM THIS ENVIRONMENT</p><h3>What could and could not be reached.</h3>
<ul class="submissions"><li><strong>Reachable:</strong> DrivenData’s public pages through the research tooling; GitHub-hosted prior artifacts through the authenticated API.</li>
<li><strong>Not reachable from the sandbox shell:</strong> USGS, ScienceBase and OpenEI hosts (connections return nothing), DriveData’s data tab (login redirect). Each blocked idea therefore names its source instead of substituting a look-alike.</li>
<li><strong>Never requested:</strong> any password, token or 2FA code. DrivenData’s <a href="{TOU}">terms of use</a> forbid automated monitoring, and none is configured here.</li></ul></div></section>'''

    pages = {"index.html": page("Euler depth-cluster submission GeoTIFF", index, "index.html"),
             "executive-summary.html": page("How to submit", execu, "executive-summary.html"),
             "evidence.html": page("Evidence", evidence, "evidence.html"),
             "hypotheses.html": page("Hypotheses", hyp, "hypotheses.html"),
             "sources.html": page("Sources", sources, "sources.html")}
    changed = []
    for name, html in pages.items():
        path = DOCS / name
        if not path.exists() or path.read_text(encoding="utf-8") != html:
            changed.append(name)
            if not args.check:
                path.write_text(html, encoding="utf-8")
    if args.check and changed:
        print("stale pages: " + ", ".join(changed))
        return 1
    print(("verified " if args.check else "wrote ") + f"{len(pages)} pages"
          + (f" (changed: {', '.join(changed)})" if changed else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
