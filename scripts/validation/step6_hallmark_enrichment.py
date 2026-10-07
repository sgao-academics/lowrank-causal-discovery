# -*- coding: utf-8 -*-
"""P0-5 fix: replace mislabelled '16 MSigDB hallmark' enrichment with a genuine
MSigDB-hallmark test computed on DISCOVERED edges only, against BOTH a uniform
null and a degree-matched (stub-matching) null.  Writes a JSON + prints a log.
"""
import csv, json, os, sys, random, collections
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')

CS = r'.'
EDGE = os.path.join(CS, 'checkpoints', 'edge_list.csv')
HALL_CANDIDATES = [
    r'{DATA_ROOT}/cancer_application/data/validation/msigdb_hallmark.json',
    r'{DATA_ROOT}/cancer_application/data/msigdb_hallmark.json',
]

def find_hall():
    for p in HALL_CANDIDATES:
        if os.path.exists(p) and os.path.getsize(p) > 0:
            return p
    # search
    for root, dirs, files in os.walk(r'{DATA_ROOT}/cancer_application'):
        for f in files:
            if 'hallmark' in f.lower() and f.endswith('.json'):
                return os.path.join(root, f)
    return None

HP = find_hall()
print('hallmark file =', HP)
if HP is None:
    sys.exit('no hallmark json found')

rows = list(csv.DictReader(open(EDGE, encoding='utf-8', errors='ignore')))
all_genes = set()
for r in rows:
    all_genes.add(r['source_gene']); all_genes.add(r['target_gene'])
D = [(r['source_gene'], r['target_gene']) for r in rows if r['source'] == 'discovered']
print('edge rows=%d  discovered=%d  genes=%d' % (len(rows), len(D), len(all_genes)))
G = sorted(all_genes)
gidx = {g: i for i, g in enumerate(G)}

# degree sequences
out_deg = collections.Counter(a for a, b in D)
in_deg = collections.Counter(b for a, b in D)

raw = json.load(open(HP, encoding='utf-8'))
sets = {}
for name, meta in raw.items():
    if isinstance(meta, dict):
        genes = meta.get('geneSymbols') or meta.get('genes') or []
    elif isinstance(meta, list):
        genes = meta
    else:
        genes = []
    s = set(genes) & all_genes
    if len(s) >= 5:
        sets[name] = s
print('hallmark sets (>=5 genes in network) =', len(sets))

n = len(D)
Ng = len(G)

# observed within-set edges (discovered only)
obs = {}
for name, s in sets.items():
    k = sum(1 for a, b in D if a in s and b in s)
    obs[name] = k

# ---- Null A: uniform random gene pairing (label permutation) -------------
rng = np.random.default_rng(42)
labels = np.array(G)
sidx = {name: np.array([gidx[g] for g in s]) for name, s in sets.items()}
obs_by_idx = {}
for name, s in sets.items():
    obs_by_idx[name] = obs[name]

# ---- Null B: degree-matched stub matching -------------------------------
# build stub lists, shuffle, pair; repeat R times
out_stubs = []
in_stubs = []
for g in G:
    out_stubs += [gidx[g]] * out_deg.get(g, 0)
    in_stubs += [gidx[g]] * in_deg.get(g, 0)
assert len(out_stubs) == len(in_stubs) == n, (len(out_stubs), n)

R = 300
nullB = {name: np.zeros(R, dtype=int) for name in sets}
nullA = {name: np.zeros(R, dtype=int) for name in sets}
oa = np.array(out_stubs); ia = np.array(in_stubs)
for t in range(R):
    # A: permute labels uniformly (keeps degree of each *position*, destroys gene identity)
    pa = rng.permutation(oa); pb = rng.permutation(ia)
    # B: stub matching - shuffle the partner lists independently (preserves each gene's degree)
    sb = rng.permutation(ia)
    for name, s in sets.items():
        sa = sidx[name]
        mask_sa = np.isin(pa, sa); mask_sb = np.isin(sb, sa)
        nullB[name][t] = int((mask_sa & mask_sb).sum())
        mask_sb2 = np.isin(pb, sa)
        nullA[name][t] = int((mask_sa & mask_sb2).sum())

out = []
for name in sets:
    K = len(sets[name])
    o = obs[name]
    eA = float(nullA[name].mean()); eB = float(nullB[name].mean())
    sdA = float(nullA[name].std()); sdB = float(nullB[name].std())
    zA = (o - eA) / sdA if sdA > 0 else 0.0
    zB = (o - eB) / sdB if sdB > 0 else 0.0
    pA = (1 + int((nullA[name] >= o).sum())) / (R + 1)
    pB = (1 + int((nullB[name] >= o).sum())) / (R + 1)
    out.append(dict(name=name, K=K, obs=o,
                    exp_uniform=eA, fold_uniform=(o / eA if eA > 0 else float('nan')), z_uniform=zA, p_uniform=pA,
                    exp_deg=eB, fold_deg=(o / eB if eB > 0 else float('nan')), z_deg=zB, p_deg=pB))

out.sort(key=lambda d: d['p_uniform'])
print('\n=== MSigDB hallmark enrichment of DISCOVERED edges (n=%d, genes=%d) ===' % (n, Ng))
print('%-34s %4s %5s %9s %8s %7s %9s %8s %7s' % ('hallmark', 'K', 'obs', 'expUnif', 'foldU', 'pU', 'expDeg', 'foldD', 'pD'))
for d in out[:20]:
    print('%-34s %4d %5d %9.2f %8.1f %7.4f %9.2f %8.1f %7.4f' % (
        d['name'][:34], d['K'], d['obs'], d['exp_uniform'], d['fold_uniform'], d['p_uniform'],
        d['exp_deg'], d['fold_deg'], d['p_deg']))
sigU = [d for d in out if d['p_uniform'] < 0.05]
sigD = [d for d in out if d['p_deg'] < 0.05]
print('\nsets with p_uniform<0.05: %d / %d' % (len(sigU), len(out)))
print('sets with p_deg<0.05    : %d / %d' % (len(sigD), len(out)))
print('max fold (uniform) = %.1f (%s)' % (max(d['fold_uniform'] for d in out), max(out, key=lambda d: d['fold_uniform'])['name']))
print('max fold (deg)     = %.1f (%s)' % (max(d['fold_deg'] for d in out), max(out, key=lambda d: d['fold_deg'])['name']))
nk = np.array([d['K'] for d in out]); nfo = np.array([d['fold_uniform'] for d in out])
print('corr(K, fold_uniform) = %.3f' % np.corrcoef(nk, nfo)[0, 1])

json.dump(out, open(os.path.join(CS, 'checkpoints', 'hallmark_enrichment_discovered.json'), 'w', encoding='utf-8'), indent=1)
print('\nwrote checkpoints/hallmark_enrichment_discovered.json')
