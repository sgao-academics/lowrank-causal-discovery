# -*- coding: utf-8 -*-
"""Give the network-gene-set enrichment claim proper statistics (audit P2-4):
up-regulation rate among network genes vs the rest, with OR + 95% CI + Fisher p.
Also verify the paper's printed 29.3/23.9 and 25.7/18.0 figures."""
import os, sys, json, csv
import numpy as np
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')
GW = r'{DATA_ROOT}/cancer_application/data/indian_oral/indian_genomewide_de.json'
EDGE = r'./checkpoints/edge_list.csv'
out = []
gw = json.load(open(GW, encoding='utf-8'))
P = gw['p']; D = gw['delta']
net = set()
for r in csv.DictReader(open(EDGE, encoding='utf-8', errors='ignore')):
    net.add(r['source_gene']); net.add(r['target_gene'])

for acc in ('GSE85195', 'GSE23558'):
    p = P[acc]; d = D[acc]
    allg = [g for g in p if g in d]
    n_net = [g for g in allg if g in net]
    n_rest = [g for g in allg if g not in net]
    a = sum(1 for g in n_net if d[g] > 0 and p[g] < 0.05)
    b = len(n_net) - a
    c = sum(1 for g in n_rest if d[g] > 0 and p[g] < 0.05)
    dd = len(n_rest) - c
    orr = (a * dd) / (b * c)
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / dd)
    tab = np.array([[a, b], [c, dd]])
    chi2, pchi, dof, exp = stats.chi2_contingency(tab, correction=False)
    fp = stats.fisher_exact(tab)[1]
    out.append('=== %s ===' % acc)
    out.append('  network genes n=%d  up&p<0.05 = %d (%.1f%%)' % (len(n_net), a, 100 * a / len(n_net)))
    out.append('  rest       n=%d  up&p<0.05 = %d (%.1f%%)' % (len(n_rest), c, 100 * c / len(n_rest)))
    out.append('  ratio = %.3f   OR = %.3f (95%% CI %.3f-%.3f)  chi2 p=%.3e  Fisher p=%.3e'
               % ((a / len(n_net)) / (c / len(n_rest)), orr, orr * np.exp(-1.96 * se), orr * np.exp(1.96 * se), pchi, fp))

# top out-degree / in-degree subsets, binomial test vs baseline
import collections
od = collections.Counter(); ideg = collections.Counter()
for r in csv.DictReader(open(EDGE, encoding='utf-8', errors='ignore')):
    if r['source'] == 'discovered':
        od[r['source_gene']] += 1; ideg[r['target_gene']] += 1
p1 = P['GSE85195']; p2 = P['GSE23558']; d1 = D['GSE85195']; d2 = D['GSE23558']
meas = [g for g in net if g in p1 and g in p2]
def bh(pv, gs):
    a = np.array([pv[x] for x in gs]); o = np.argsort(a); m = len(a)
    q = np.empty(m); prev = 1.0
    for i in range(m - 1, -1, -1):
        idx = o[i]; prev = min(prev, a[idx] * m / (i + 1)); q[idx] = prev
    return {gs[i]: q[i] for i in range(m)}
q1 = bh(p1, meas); q2 = bh(p2, meas)
both = set(g for g in meas if q1[g] < 0.05 and q2[g] < 0.05 and d1[g] > 0 and d2[g] > 0)
out.append('\nboth-cohort up (FDR5) set = %d ; baseline %d/%d = %.1f%%' % (len(both), len(both), len(meas), 100 * len(both) / len(meas)))
base = len(both) / len(meas)
for nm, cnt in (('out-degree', od), ('in-degree', ideg)):
    order = sorted(meas, key=lambda g: -cnt.get(g, 0))
    for k in (50, 100, 200):
        hit = sum(1 for g in order[:k] if g in both)
        pb = stats.binomtest(hit, k, base, alternative='greater').pvalue
        out.append('  top-%3d %s: %3d/%3d = %.1f%%  binom p=%.3e (baseline %.1f%%)'
                   % (k, nm, hit, k, 100 * hit / k, pb, 100 * base))
open(r'./_enrich_test.txt', 'w', encoding='utf-8').write('\n'.join(out))
print('\n'.join(out))
