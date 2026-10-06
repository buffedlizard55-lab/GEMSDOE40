"""Harvest prior scored rasters (public sibling repos) for education/audit only.

Every entry records the owner-reported live score attached to that file by the sibling
project's own site/ledger.  These files are NEVER re-submitted; they are used to
validate the local proxy instrument and to audit uniqueness.
"""
import json, subprocess, sys, os
from concurrent.futures import ThreadPoolExecutor

ITEMS = [
 # (repo, path, local_name, live_score)
 ("GEMSDOE24","inputs/calibration/13gems_20261001_r13-lattice-s5_v2_nan-outside.tif","p13-lattice-s5-v2.tif",0.0904),
 ("GEMSDOE24","inputs/calibration/8GEMSDOE_Hedge-v2_submission.tif","p08-hedge-v2.tif",0.1563),
 ("GEMSDOE24","inputs/calibration/gemsdoe-ens12-adopted-7f00890a.tif","p01-ens12-7f00890a.tif",0.1563),
 ("GEMSDOE24","inputs/calibration/gemsdoe9-PLACEHOLDER-2314b599.tif","p09-2314b599.tif",0.0107),
 ("GEMSDOE24","inputs/calibration/gems10-h25-ctx-ridge-20260927T232947704150Z-6452ae1d00.tif","p10-h25-ctx-ridge.tif",0.1280),
 ("GEMSDOE24","inputs/calibration/gems10-h28-dotted-ridge-20260928T020256236880Z-6452ae1d00.tif","p10-h28-dotted-ridge.tif",0.1839),
 ("GEMSDOE22","docs/downloads/gems22-h23-a-dti-optimal-emission-6pct-20261002-e2ec4b49-allfinite.tif","p22-h23a-6pct.tif",0.1002),
 ("GEMSDOE22","docs/downloads/gems22-h23-b-dti-optimal-emission-10pct-20261002-86176698-allfinite.tif","p22-h23b-10pct.tif",0.0748),
 ("GEMSDOE22","assets/lb_anchors/lb_anchor_h19-5_0.1922.tif","p19-h19-5-1922.tif",0.1922),
 ("GEMSDOE22","assets/lb_anchors/lb_anchor_h19-4_0.1894.tif","p19-h19-4-1894.tif",0.1894),
 ("GEMSDOE22","assets/lb_anchors/lb_anchor_h16-1_0.1855.tif","p16-h16-1-1855.tif",0.1855),
 ("GEMSDOE24","docs/downloads/gems24-h25-1-dotted-h19-5-d1-5-20261002-989f59505db1-allfinite.tif","p24-d15-2477.tif",0.2477),
 ("GEMSDOE24","docs/downloads/gems24-h25-1-dotted-h19-5-d2-8-20261002-e56ea318af89-allfinite.tif","p24-d28-2600.tif",0.2600),
 ("GEMSDOE24","docs/downloads/gems24-reference-h19-5-20261002-80d47e1ab2ee-allfinite.tif","p19-h19-5-solid.tif",0.1922),
 ("GEMSDOE27","docs/downloads/gems27-topo-gap-closure-t-v2-on-d1-5-20261002-5512495c6bd1-allfinite.tif","p27-topogap-d15.tif",0.2449),
 ("GEMSDOE28","docs/downloads/gems28-h27-4-r1-solo-d2-8-20261003-8acb75e1f2cc-allfinite.tif","p28-h27-4-d28.tif",0.2708),
 ("GEMSDOE30","docs/downloads/gemsdoe30-d28-poisson300m-offcat-44090-20261003T233156Z-91eae1ca.tif","p30-d28-poisson300-offcat.tif",0.2600),
 ("GEMSDOE30","docs/downloads/gemsdoe30-sgmc-hedge-d10-85k-20261003-ac08b41e-zeros.tif","p30-sgmc-hedge-d10-85k.tif",None),
 ("GEMSDOE33","archive/legacy_candidates/gems33-c0-scored-reference-20261004-89bf5b9a2fea.tif","p33-c0-scored-ref.tif",None),
 ("GEMSDOE34","docs/downloads/h34-scatter-q50-arr-matched-20261004T223317Z.tif","p34-scatter-q50.tif",0.0778),
 ("GEMSDOE37","docs/downloads/gemsdoe37-h6-physics-dotted-80k-20261005T055000Z-0bef9211631c.tif","p37-h6-phys-dotted-80k.tif",None),
 ("GEMSDOE39","docs/downloads/gemsdoe39-h40-e-disc-h40e-30k-zeros.tif","p39-h40e-30k.tif",None),
]

def fetch(item):
    repo, path, name, score = item
    dest = f"ref/prior/{name}"
    if os.path.exists(dest) and os.path.getsize(dest) > 200000:
        return name, "cached", os.path.getsize(dest)
    r = subprocess.run(["gh","api",f"/repos/buffedlizard55-lab/{repo}/contents/{path}",
                        "-H","Accept: application/vnd.github.raw"], capture_output=True, timeout=600)
    if r.returncode != 0 or len(r.stdout) < 100000:
        return name, "FAILED", len(r.stdout)
    open(dest,"wb").write(r.stdout)
    return name, "ok", len(r.stdout)

with ThreadPoolExecutor(max_workers=5) as ex:
    for name, status, size in ex.map(fetch, ITEMS):
        print(f"{status:8s} {size:10d} {name}")

json.dump({n: {"repo": r, "path": p, "live_score": s} for r,p,n,s in ITEMS},
          open("ref/meta/prior_index.json","w"), indent=1)
