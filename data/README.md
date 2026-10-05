# Competition rasters

| file | sha256 | source |
|---|---|---|
| `sample_submission.tif` | `2176d08e485aa2cd2860ce8df539db4faf4d76163b38a4dd8c30a40454d35cbc` | DrivenData, bridged via GEMSDOE10 |
| `labels.tif` / `existing_faults.tif` | `7ba308ccdc4418b31a178f4f1ef21aaa6e152e4028f2f6f64b01f7eb25ae4093` | DrivenData |
| `training_features.tif` (gitignored, 419 MB) | `4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5` | DrivenData; assemble with `bash scripts/download_competition_data.sh` from `data/bridge/` parts or download from https://www.drivendata.org/competitions/306/competition-doe-gems/data/ (login) |

Official data tab is login-walled. This sandbox used the sha256-pinned 5GEMSDOE data-bridge.
