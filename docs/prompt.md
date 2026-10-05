# Standing brief for this repository (read before every session)

> Transcribed from the owner's brief. Repeated paragraphs that appear several times in the original
> message are consolidated here **once**; all substantive requirements are preserved, and the
> highest-urgency block is verbatim.

## Highest urgency (verbatim)

> THE FOLLOWING IS THE HIGHEST URGENCY AND MUST BE FOLLOWED!
>
> MUST GENERATE A UNIQUE TIF SUBMISSION FOR THE COMPETITION. DO NOT COPY A PREVIOUS SUBMISSION
> UNLESS IT'S FOR LEARNING AND EDUCATION. BUT WE MUST GENERATE A UNIQUE TIF SUBMISSION.
>
> There should be an easy to download submission tif file as described by the prompt. Read the
> entire prompt.
>
> Euler deconvolution depth-clustering. Build this from potential-field depth estimation, not
> gradient thresholding: run Euler deconvolution (Reid, Allsop, Granser, Millett, and Somerton,
> Geophysics, 1990) across the magnetic and gravity layers with the structural index for a
> fault-like contact, producing a cloud of depth-labeled solution points rather than a single edge
> map. Convert that cloud into a continuous raster via kernel-density estimation of solution density
> per pixel, weighted so tight clusters of shallow, mutually-consistent solutions score higher than
> scattered or deep ones, since a real near-surface fault produces the former and noise produces the
> latter. Normalize to [0,1] and write a single-band float32 GeoTIFF in EPSG:32611, 100 m
> resolution, matching the sample submission's exact shape and geotransform, NaN only outside the
> valid footprint — then, before presenting it for download, hash and correlate it against every
> prior submission's raw output and refuse to call it new if it's a near-duplicate, since
> depth-clustering should produce a visibly different spatial pattern than any gradient or curvature
> candidate already made.

## Competition and data

* Competition: **The Geologic Enhanced Mapping System (GEMS) Prize Challenge** — DrivenData
  competition 306, U.S. DOE Office of Geothermal, <https://www.drivendata.org/competitions/306/competition-doe-gems/>.
* Problem description: <https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/>;
  additional resources: <https://www.drivendata.org/competitions/306/competition-doe-gems/page/968/>.
* Data tab (login-walled): <https://www.drivendata.org/competitions/306/competition-doe-gems/data/>.
* Rules PDF: <https://docs.nlr.gov/docs/fy26osti/96647.pdf>.
* Leaderboard: <https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/>;
  current #1 = 0.3195; this programme's best owner-reported score = 0.2778.
* Reference solution: <https://github.com/drivendataorg/gems-prize-reference-solution>.
* Submission contract: a single-band GeoTIFF (or a .zip containing one), 100 m, EPSG:32611, same
  shape and geotransform as the sample submission, values in [0, 1] inside the valid footprint.
  A short Note field is required when uploading (unique name + short comment).

## How work must be done

* **Generate a unique TIF submission** every session; never re-upload a previous submission's bytes.
* **Hash and correlate** every artifact against all prior submissions' raw outputs; refuse to call
  it new if it is a near-duplicate.
* **Before implementing**, generate 3–5 candidate geological hypotheses not yet tried, each naming:
  the specific layer(s) involved, the physical signature targeted, why it should catch a fault
  missing from the USGS/INGENIOUS catalogue rather than one already in it, and how it differs from
  anything already implemented in the repository. Rank them by expected DTI improvement and
  implementation cost.
* **Validate the top candidate on the spatially blocked holdout set before touching a weekly
  submission slot.** Do not spend a submission slot on an idea that has not beaten the current
  holdout best.
* If a candidate **cannot be validated without new external data**, name the specific free, official
  source needed and check that it is obtainable **before** proposing the idea as viable.
* Work **line by line**, verifying against official trusted sources, **with links for manual
  review**. **No manual input** — work autonomously. **Flag irregularities** for review.
  **No hallucinations.**
* The site must make the submission file **obvious and one click away** from the landing page, and
  must include an **executive-summary subpage** explaining exactly how to make a submission —
  including the fix for the "Predicted values must be in range [0, 1]" error seen on the portal.
* Keep the **Core Values** central: *Maximize P(Win)* and *Own the Outcome*.
* Maintain a full, auditable list (data sources, hypotheses, measurements, next steps) with official
  verified links so that nothing has to be manually re-checked, and keep it up to date.

## Standing results context (owner-reported)

The programme's earlier repositories (GEMSDOE … GEMSDOE39) produced the score history summarised on
the [leaderboard analysis page](leaderboard-analysis.html). Best owner-reported score: **0.2778**
(GEMSDOE32 `h33-h33-2-b2`). The best score currently on the competition leaderboard is **0.3195**.

## Deliverables of this repository

1. A one-click submission GeoTIFF in `docs/downloads/` (served by the GitHub Pages site).
2. A validated-verdict banner: what to upload and what *not* to upload, with the measurements behind
   both statements.
3. The auditable evidence trail: band inventory, submission-format verification, instrument
   verification, candidate sweeps, novelty/hash audit.
4. Named, ranked next hypotheses — including the ones that were refuted, and why.
