"""sha256 pins for every external byte this repository reads.

A pin proves the bytes are the ones that were inspected; it does **not** prove organizer
authentication (the organizer's data tab is login-walled).  The pins were established by reading
the files from byte-identical GitHub mirrors and are re-checked by ``scripts/inspect_data.py``.
"""

PINNED_FILES: dict[str, dict] = {
    "training_features.tif": dict(
        bytes=418_912_844,
        sha256="4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5",
        canonical_name="training_features.tif",
        mirror="github.com/buffedlizard55-lab/GEMSDOE data/bridge/gems-geodawn-numerical-features.tif.part-000..004",
        note="19 float32 bands; band tags carry the organiser's own band names"),
    "existing_faults.tif": dict(
        bytes=425_830,
        sha256="7ba308ccdc4418b31a178f4f1ef21aaa6e152e4028f2f6f64b01f7eb25ae4093",
        canonical_name="labels.tif",
        mirror="github.com/buffedlizard55-lab/GEMSDOE24 data/bridge/labels.tif",
        note="USGS Quaternary + INGENIOUS catalogue, int8, -1 outside the footprint"),
    "example_submission.tif": dict(
        bytes=1_599_597,
        sha256="2176d08e485aa2cd2860ce8df539db4faf4d76163b38a4dd8c30a40454d35cbc",
        canonical_name="sample_submission.tif",
        mirror="github.com/buffedlizard55-lab/GEMSDOE24 data/bridge/sample_submission.tif",
        note="geometry template; 5,167,373-pixel footprint, 0/1 inside, NaN outside"),
    "derived_sgmc_faults_100m_u8.tif": dict(
        bytes=198_602,
        sha256="643cbe992ef4ba37588fb469163ed8291e3ceb23d6c1f78a3cfaa462430c2da0",
        canonical_name="derived_sgmc_faults_100m_u8.tif",
        mirror="github.com/buffedlizard55-lab/GEMSDOE24 data/external/derived_sgmc_faults_100m_u8.tif",
        note="current owner-derived SGMC validation proxy only; not organizer truth and distinct from the archived H2-B proxy bytes"),
}
