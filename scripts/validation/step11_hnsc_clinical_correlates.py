# -*- coding: utf-8 -*-
"""External validation #5: is the core uniform across head-and-neck sub-sites,
and does it say anything about extent or outcome?

Everything here is one population on one platform (TCGA-HNSC), so the only
comparison that carries information is between sub-sites, and that comparison
has to be made against a matched null: the oral cavity and the larynx differ
transcriptome-wide, so a random gene set of the same size already separates
them (it separates them downwards).  Section S14 reports the contrast, the
matched null, the within-tissue-source-site control, the per-gene values, and
the null results for pathological stage and overall survival.

Data (all public, not redistributed with the repository):
  expression  {DATA_ROOT}/cancer_application/data/TCGA_HNSC_HiSeqV2.tsv
              (UCSC Xena, TCGA.HNSC.sampleMap/HiSeqV2)
  clinical    {DATA_ROOT}/cancer_application/data/TCGA_HNSC_clinical.tsv
              (UCSC Xena, TCGA.HNSC.sampleMap/HNSC_clinicalMatrix --
               anatomic_neoplasm_subdivision, pathologic_stage,
               hpv_status_by_p16_testing)
  survival    {DATA_ROOT}/cancer_application/data/validation/pancan_os/HNSC_os.json
              (UCSC Xena pan-cancer survival release)
"""
import csv
import json
import os
import sys

import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
T = r'{DATA_ROOT}/cancer_application/data/TCGA_HNSC_HiSeqV2.tsv'
CLIN = r'{DATA_ROOT}/cancer_application/data/TCGA_HNSC_clinical.tsv'
OSF = r'{DATA_ROOT}/cancer_application/data/validation/pancan_os/HNSC_os.json'
EDGE = r'./checkpoints/edge_list.csv'

NINE = ['STAT1', 'FOSL1', 'CXCL8', 'MMP9', 'E2F1', 'CXCL10', 'JUN', 'RELA',
        'ETS1']
ORAL = {'Oral Cavity', 'Buccal Mucosa', 'Alveolar Ridge', 'Floor of mouth',
        'Oral Tongue', 'Hard Palate', 'Lip'}
LARYNX = {'Larynx', 'Hypopharynx'}
out = []

# ------------------------------------------------------------------ expression
hdr = open(T, encoding='utf-8', errors='ignore').readline().rstrip('\n').split('\t')
samples = hdr[1:]
keep = [i for i, s in enumerate(samples) if s[13:15] == '01']
kept = [samples[i] for i in keep]
tss = np.array([s[5:7] for s in kept])
pid = [s[:12] for s in kept]
names, rows = [], []
with open(T, encoding='utf-8', errors='ignore') as f:
    f.readline()
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) == len(samples) + 1:
            names.append(p[0])
            rows.append(np.array(p[1:], dtype=float)[keep])
M = np.array(rows)
M = (M - M.mean(1, keepdims=True)) / (M.std(1, keepdims=True) + 1e-12)
gi = {g: i for i, g in enumerate(names)}
out.append('TCGA-HNSC tumours=%d, genes=%d' % (len(kept), len(names)))
core = [g for g in NINE if g in gi]
out.append('core genes measurable: %d of 9  (%s absent)'
           % (len(core), ','.join(g for g in NINE if g not in gi)))
K = len(core)
nine = M[[gi[g] for g in core]].mean(0)

# -------------------------------------------------------------------- clinical
clin = {}
for r in csv.DictReader(open(CLIN, encoding='utf-8', errors='replace'),
                        delimiter='\t'):
    clin[r['sampleID'][:12]] = r
site = np.array([clin[p]['anatomic_neoplasm_subdivision'] if p in clin else ''
                 for p in pid])
hpv = np.array([clin[p]['hpv_status_by_p16_testing'] if p in clin else ''
                for p in pid])
male = np.array([clin[p]['gender'] == 'MALE' if p in clin else False
                 for p in pid])
o, l = np.isin(site, list(ORAL)), np.isin(site, list(LARYNX))


def cmp(a, b, tag):
    p = stats.mannwhitneyu(a, b, alternative='two-sided')[1]
    d = (a.mean() - b.mean()) / np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    out.append('  %-40s n=%3d vs %3d  %+.4f vs %+.4f  p=%.3e  d=%+.2f'
               % (tag, len(a), len(b), a.mean(), b.mean(), p, d))
    return p


out.append('\n=== sub-site contrast (eight-gene core score) ===')
out.append('  oral-cavity n=%d (%.1f%% of tumours), larynx+hypopharynx n=%d'
           % (o.sum(), 100 * o.mean(), l.sum()))
cmp(nine[o], nine[l], 'oral vs larynx, all tumours')
cmp(nine[o & (hpv == 'Negative')], nine[l & (hpv == 'Negative')],
    'p16-negative only')
cmp(nine[o & male], nine[l & male], 'men only')
out.append('  sub-site means:')
for k in sorted(set(site), key=lambda k: -int((site == k).sum())):
    if k:
        out.append('     %-18s n=%3d  %+.3f' % (k, (site == k).sum(),
                                                nine[site == k].mean()))

# ------------------------------------------------- matched null, TSS control
mu = {t: M[:, tss == t].mean(1, keepdims=True) for t in set(tss)}
Mc = M.copy()
for t in set(tss):
    Mc[:, tss == t] -= mu[t][:, 0][:, None]
rng = np.random.default_rng(20261008)
pool = np.array([i for i, g in enumerate(names) if g not in set(core)])
out.append('\n=== matched null: 4000 random gene sets of the same size ===')
for tag, mat in (('raw', M), ('within-TSS-centred', Mc)):
    v = mat[[gi[g] for g in core]].mean(0)
    obs = v[o].mean() - v[l].mean()
    null = np.empty(4000)
    for b in range(4000):
        gs = pool[rng.choice(len(pool), K, replace=False)]
        w = mat[gs].mean(0)
        null[b] = w[o].mean() - w[l].mean()
    out.append('  [%s] observed %+.4f | null mean %+.4f sd %.4f range [%+.4f,'
               ' %+.4f] | draws >= observed: %d | z=%+.2f'
               % (tag, obs, null.mean(), null.std(), null.min(), null.max(),
                  int((null >= obs).sum()),
                  (obs - null.mean()) / null.std()))

out.append('\n=== per gene, oral minus larynx (z units) ===')
for g in NINE:
    if g in gi:
        v = M[gi[g]]
        out.append('  %-7s %+.4f  p=%.3e' % (g, v[o].mean() - v[l].mean(),
                                             stats.mannwhitneyu(v[o], v[l])[1]))

# ------------------------------------------------------------ extent, outcome
ORD = {'T1': 1, 'T2': 2, 'T3': 3, 'T4': 4}


def lvl_T(r):
    return ORD.get(r['pathologic_T'].rstrip('abcd'))


def lvl_N(r):
    return {'N0': 0, 'N1': 1, 'N2': 2, 'N3': 3}.get(r['pathologic_N'][:2])


def lvl_S(r):
    for k, n in (('Stage IV', 4), ('Stage III', 3), ('Stage II', 2),
                 ('Stage I', 1)):
        if r['pathologic_stage'].startswith(k):
            return n


out.append('\n=== disease extent (Kruskal-Wallis on the score) ===')
for name, fn in (('pathologic T', lvl_T), ('pathologic N', lvl_N),
                 ('overall stage', lvl_S)):
    gr = {}
    for i, p in enumerate(pid):
        if p in clin and fn(clin[p]) is not None:
            gr.setdefault(fn(clin[p]), []).append(nine[i])
    out.append('  %-16s p=%.4f  %s' % (name, stats.kruskal(
        *[np.array(gr[k]) for k in sorted(gr)])[1], '  '.join(
            '%s:%+.2f(n=%d)' % (k, np.mean(gr[k]), len(gr[k]))
            for k in sorted(gr))))

OS = {x['patient']: x for x in json.load(open(OSF, encoding='utf-8'))}


def logrank_k(t, e, g):
    """k-group log-rank, chi-square with k-1 degrees of freedom."""
    t, e, g = np.asarray(t, float), np.asarray(e, bool), np.asarray(g)
    ks = sorted(set(g.tolist()))
    O = {k: 0.0 for k in ks}
    Ex = {k: 0.0 for k in ks}
    V = np.zeros((len(ks), len(ks)))
    for ti in np.unique(t[e]):
        R = t >= ti
        n = np.array([float((R & (g == k)).sum()) for k in ks])
        N = n.sum()
        if N < 2:
            continue
        d = float((e & (t == ti)).sum())
        dk = np.array([float((e & (t == ti) & (g == k)).sum()) for k in ks])
        pr = n / N
        for a, k in enumerate(ks):
            O[k] += dk[a]
            Ex[k] += d * pr[a]
        if N - d > 0:
            V += d * (np.diag(pr) - np.outer(pr, pr)) * (N - d) / (N - 1)
    z = np.array([O[k] - Ex[k] for k in ks])
    Vr = V[:-1, :-1]
    chi2 = float(z[:-1] @ np.linalg.pinv(Vr) @ z[:-1]) if Vr.size else 0.0
    return O, Ex, chi2, float(stats.chi2.sf(chi2, max(len(ks) - 1, 1)))


def logrank2(t, e, hi):
    t, e, hi = np.asarray(t, float), np.asarray(e, bool), np.asarray(hi, bool)
    O1 = E1 = V1 = 0.0
    for ti in np.unique(t[e]):
        R = t >= ti
        n1, n2 = int((hi & R).sum()), int((~hi & R).sum())
        N = n1 + n2
        if N < 2:
            continue
        d = int((e & (t == ti)).sum())
        d1 = int((e & (t == ti) & hi).sum())
        pr = n1 / N
        O1 += d1
        E1 += d * pr
        V1 += d * pr * (1 - pr) * (N - d) / (N - 1)
    chi2 = (O1 - E1) ** 2 / V1 if V1 > 0 else 0.0
    return O1, E1, float(stats.chi2.sf(chi2, 1))


sv = [(i, OS[p]['os_months'], OS[p]['os_status'].startswith('1'))
      for i, p in enumerate(pid) if p in OS and OS[p]['os_months'] > 0]
idx = np.array([i for i, _, _ in sv])
TT = np.array([m for _, m, _ in sv], float)
DD = np.array([e for _, _, e in sv], bool)
out.append('\n=== overall survival (n=%d tumours, %d deaths) ==='
           % (len(TT), int(DD.sum())))
out.append('  median split  p=%.3f   oral-cavity only p=%.3f'
           % (logrank2(TT, DD, (nine[idx] > np.median(nine[idx])))[2],
              logrank2(TT[o[idx]], DD[o[idx]],
                       (nine[idx][o[idx]] > np.median(nine[idx][o[idx]])))[2]))
lo = np.quantile(nine[idx], 1 / 3)
hi_q = np.quantile(nine[idx], 2 / 3)
ter = np.digitize(nine[idx], [lo, hi_q])
out.append('  tertiles      p=%.3f' % logrank_k(TT, DD, ter)[3])
for tag, sel in (('p16-positive', hpv[idx] == 'Positive'),
                 ('p16-negative', hpv[idx] == 'Negative')):
    if sel.sum() > 15:
        out.append('  %-13s n=%3d  p=%.3f'
                   % (tag, sel.sum(),
                      logrank2(TT[sel], DD[sel],
                               (nine[idx][sel] > np.median(nine[idx][sel])))[2]))

# --------------------------------------------------------------- nine genes
edge = list(csv.DictReader(open(EDGE, encoding='utf-8')))
outd, ind = {}, {}
for r in edge:
    if r['source'] == 'discovered':
        s, t = r['source_gene'].upper(), r['target_gene'].upper()
        outd[s] = outd.get(s, 0) + 1
        ind[t] = ind.get(t, 0) + 1
out.append('\n=== degrees of the nine in the released estimator graph ===')
for g in NINE:
    out.append('  %-7s out=%4d  in=%4d' % (g, outd.get(g, 0), ind.get(g, 0)))

print('\n'.join(out))
