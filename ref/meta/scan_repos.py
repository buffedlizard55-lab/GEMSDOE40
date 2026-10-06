import json, subprocess, sys, os
repos = [l.split('\t')[0] for l in open('ref/meta/repos.tsv') if l.strip().startswith(('GEMSDOE',))]
repos = [r for r in repos if r not in ('GEMSDOE43','GEMSDOE44','GEMSDOE45','GEMSDOE46','GEMSDOE47','GEMSDOE48','GEMSDOE49','GEMSDOE50','GEMSDOE51','GEMSDOE52','GEMSDOE53','GEMSDOE54')]
out = {}
for r in repos:
    try:
        t = subprocess.run(['gh','api',f'/repos/buffedlizard55-lab/{r}/git/trees/HEAD?recursive=1'],
                           capture_output=True, text=True, timeout=120).stdout
        d = json.loads(t)
        blobs = [(e['path'], e['size']) for e in d.get('tree',[]) if e['type']=='blob']
        out[r] = blobs
    except Exception as e:
        print('ERR', r, e, file=sys.stderr)
json.dump(out, open('ref/meta/repo_trees.json','w'))
print('repos scanned', len(out))
tifs = [(r,p,s) for r,bs in out.items() for p,s in bs if p.lower().endswith('.tif')]
print('total .tif blobs', len(tifs))
from collections import Counter
print(Counter(r for r,_,_ in tifs).most_common(30))
