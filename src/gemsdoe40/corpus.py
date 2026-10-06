"""The prior-submission corpus used for calibration and novelty auditing.

Every entry is a *public* GeoTIFF published by a sibling GEMSDOE repository of
the same owner, with the score that the owner/user reported for it.  The scores
are owner-reported, not organiser receipts: the public leaderboard does not
expose file names or hashes, so the (file, score) pairing is an assumption that
is disclosed wherever these numbers are used.

The corpus is fetched on demand (``scripts/fetch_prior_corpus.sh``) into
``data/prior/`` from the pinned repository paths below; nothing here is
committed.  ``sha256`` is of the *file as published in that repository*.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PRIOR_DIR = ROOT / "data" / "prior"


@dataclass(frozen=True)
class Prior:
    key: str
    score: float
    repo: str
    path: str
    sha256: str
    emitted_px: int | None = None
    notes: str = ""
    #: False when the artifact was constructed from the surrogate catalogue
    #: itself (SGMC), which makes surrogate scoring circular.
    use_for_calibration: bool = True


#: Emitted pixel counts are measured from the files by ``scripts/study_priors.py``;
#: the values recorded here are the owner-reported ones where available.
PRIORS: tuple[Prior, ...] = (
    Prior("h33-2-b2", 0.2778, "GEMSDOE32",
          "docs/downloads/gemsdoe32-h33-h33-2-b2-20261004T220000Z-e5eb6e7e-nan.tif",
          "baeae3219bba6a19bc8cfe79224e28555c50184b7c1ca65c2a822f54807c77dd", 37654,
          "0.2708 base minus every dot within 2 px of the catalogue; site says UNSCORED"),
    Prior("h27-4-r1-solo-d2-8", 0.2708, "GEMSDOE28",
          "docs/downloads/gems28-h27-4-r1-solo-d2-8-20261003-8acb75e1f2cc-nan.tif", "", 40199,
          "dot lattice, min separation 2*sqrt(2) px"),
    Prior("h27-4-solo-d28-mirror", 0.2708, "GEMSDOE31",
          "docs/downloads/gemsdoe31-h27-4-solo-d28-20261004-8acb75e1-nan.tif", "", 40199,
          "mirror of the same candidate; duplicate for calibration"),
    Prior("h32-1-prethin-tip-euler-d2-8", 0.2649, "GEMSDOE28",
          "docs/downloads/gems28-h32-1-prethin-tip-euler-d2-8-20261003-31e35eee884e-nan.tif", "", 42294,
          "tip-protracted Euler arm added to a dotted base - closest prior art to H4"),
    Prior("h33d-analog-tip-stepover-r30", 0.2632, "GEMSDOE33",
          "docs/downloads/GEMSDOE33-h33d-analog-tip-stepover-r30-20261004-cb490425926e-nanoutside.tif",
          "", 41865, "bounded transfer + Euler corroboration + flank prune"),
    Prior("dotted-h19-5-d1-5", 0.2477, "GEMSDOE24",
          "docs/downloads/gems24-h25-1-dotted-h19-5-d1-5-20261002-989f59505db1-nan.tif", "", 60069,
          "h19-5 dot-thinned at min separation 1.5 px"),
    Prior("topo-gap-closure-v2-d1-5", 0.2449, "GEMSDOE27",
          "docs/downloads/gems27-topo-gap-closure-t-v2-on-d1-5-20261002-5512495c6bd1-nan.tif", "", 61328,
          "topographic gap closure at min separation 1.5 px"),
    Prior("dotted-h19-5-d2-8", 0.2600, "GEMSDOE25",
          "docs/downloads/gems25-dotted-h19-5-d2-8-20261002-e56ea318af89-nan.tif", "", 44090,
          "h19-5 dot-thinned at min separation 2*sqrt(2) px"),
    Prior("historical-d28-mirror", 0.2600, "GEMSDOE29",
          "docs/downloads/gemsdoe29-historical-d28-20261002-e56ea318af89-nan.tif", "", 44090,
          "byte-identical duplicate of the previous entry"),
    Prior("d28-poisson300m-44090", 0.2600, "GEMSDOE30",
          "docs/downloads/gemsdoe30-d28-poisson300m-offcat-44090-20261003T233156Z-91eae1ca-nanoutside.tif",
          "", 44090, "off-catalogue Poisson-300 m emission of the same family"),
    Prior("h19-4", 0.1894, "19GEMSDOE",
          "docs/downloads/gems19-h19-4-multiline-corroborated-openness-thermal-pop-20260930-691e4dfa-nan.tif",
          "", 123779, "multiline corroborated openness"),
    Prior("h19-5", 0.1922, "19GEMSDOE",
          "docs/downloads/gems19-h19-5-powerlaw-budget-multiline-corroborated-20260930-e27054cf-nan.tif",
          "", 121131, "powerlaw budget, min separation 1 px"),
    Prior("h16-1-topo-geophys-ridges", 0.1855, "16GEMSDOE",
          "docs/downloads/gems16-h16-1-topo-geophys-baseline-ridges-20260930-df20f65e-nan.tif", "", 123939,
          "topographic + geophysical ridge baseline"),
    Prior("h18-3a-x-complexity", 0.0976, "16GEMSDOE",
          "docs/downloads/gems16-h18-3a-topo-geophys-x-complexity-prior-20260930-c502dfab-nan.tif", "", 122120,
          "cross-complexity prior"),
    Prior("h18-4-usgs-geologic-map-gap", 0.0360, "16GEMSDOE",
          "docs/downloads/gems16-h18-4-usgs-geologic-map-faults-gap-20260930-aef8f42c-nan.tif", "", None,
          "built from SGMC state-map faults: circular against the SGMC surrogate",
          use_for_calibration=False),
    Prior("h28-dotted-ridge", 0.1839, "GEMSDOE10",
          "docs/downloads/gems10-h28-dotted-ridge-20260928T020256236880Z-6452ae1d00.tif", "", 69281,
          "dotted ridge"),
    Prior("h25-ctx-ridge", 0.1280, "GEMSDOE10",
          "docs/downloads/gems10-h25-ctx-ridge-20260927T232947704150Z-6452ae1d00.tif", "", 174232,
          "context ridge"),
    Prior("h20-dem10-scarp-thin", 0.0921, "GEMSDOE10",
          "docs/downloads/gems10-h20-dem10-scarp-thin-20260927T155223039488Z-ffc91a1686.tif", "", 153957,
          "10 m DEM scarp thinning"),
    Prior("h16-continuation", 0.0461, "GEMSDOE10",
          "docs/downloads/gems10-h16-continuation-20260927T065521077735Z-3431b83c7c.tif", "", 310042,
          "continuation arm"),
)


def local_path(p: Prior) -> Path:
    return PRIOR_DIR / p.key / Path(p.path).name


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_field(path: str | Path, footprint: np.ndarray | None = None):
    """Read a prior GeoTIFF as float, NaN-safe, outside-footprint -> 0."""
    import numpy as np
    import rasterio

    with rasterio.open(path) as s:
        a = s.read(1).astype("float64")
        nd = s.nodata
    if nd is not None:
        a = np.where(a == nd, 0.0, a)
    a = np.where(np.isfinite(a), a, 0.0)
    if footprint is not None:
        a = np.where(footprint, a, 0.0)
    return a
