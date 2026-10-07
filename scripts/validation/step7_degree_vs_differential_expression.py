# -*- coding: utf-8 -*-
"""Unbiased test (audit gap E/F): is the network's degree ranking informative about
which genes are consistently up-regulated in BOTH Indian cohorts?
Turns P0-1 into a positive, and asks whether the network adds anything over a
plain gene list."""
import os, sys, json, csv, collections
import numpy as np
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')

GW = r'{DATA_ROOT}/cancer_application/data/indian_oral/indian_genomewide_de.json'
EDGE = r'./checkpoints/edge_list.csv'
out = []

gw = json.load(open(GW, encoding='utf-8'))
P1 = gw['p']['GSE85195']; P2 = gw['p']['GSE23558']
D1 = gw['delta']['GSE85195']; D2 = gw['delta']['GSE23558']

rows = list(csv.DictReader(open(EDGE, encoding='utf-8', errors='ignore')))
out_deg = collections.Counter()
in_deg = collections.Counter()
genes = set()
for r in rows:
    genes.add(r['source_gene']); genes.add(r['target_gene'])
    if r['source'] == 'discovered':
        out_deg[r['source_gene']] += 1
        in_deg[r['target_gene']] += 1
genes = sorted(genes)
out.append('network genes = %d ; discovered edges = %d' % (len(genes), sum(out_deg.values())))

meas = [g for g in genes if g in P1 and g in P2]
out.append('measurable in both cohorts = %d' % len(meas))

def bh(pv, gs):
    a = np.array([pv[x] for x in gs]); o = np.argsort(a); m = len(a)
    q = np.empty(m); prev = 1.0
    for i in range(m - 1, -1, -1):
        idx = o[i]; prev = min(prev, a[idx] * m / (i + 1)); q[idx] = prev
    return {gs[i]: q[i] for i in range(m)}

q1 = bh(P1, meas); q2 = bh(P2, meas)
both_up = [g for g in meas if q1[g] < 0.05 and q2[g] < 0.05 and D1[g] > 0 and D2[g] > 0]
up_only = [g for g in meas if q1[g] < 0.05 and q2[g] < 0.05 and D1[g] > 0]
out.append('both-cohort FDR<0.05 & up = %d' % len(up_only))
out.append('both-cohort (any direction) = %d' % len([g for g in meas if q1[g] < 0.05 and q2[g] < 0.05]))

# Fisher combined
fisher = {}
for g in meas:
    chi = -2 * (np.log(P1[g]) + np.log(P2[g]))
    fisher[g] = stats.chi2.sf(chi, 4)

od = np.array([out_deg.get(g, 0) for g in meas], float)
neglog = np.array([-np.log10(fisher[g]) for g in meas], float)
rho, prho = stats.spearmanr(od, neglog)
out.append('\nSpearman(out-degree, -log10 Fisher p) = %.3f  p=%.2e  (n=%d)' % (rho, prho, len(meas)))

setA = set(both_up)
mA = np.array([out_deg.get(g, 0) for g in meas if g in setA], float)
mB = np.array([out_deg.get(g, 0) for g in meas if g not in setA], float)
u, pu = stats.mannwhitneyu(mA, mB, alternative='greater')
out.append('both-up genes: mean out-degree %.1f (median %.0f) vs rest %.1f (median %.0f)  MWU p=%.2e'
           % (mA.mean(), np.median(mA), mB.mean(), np.median(mB), pu))

# top-k by out-degree -> what fraction are both-up?
order = sorted(meas, key=lambda g: -out_deg.get(g, 0))
base = len(both_up) / len(meas)
out.append('baseline share both-up = %.1f%%' % (100 * base))
for k in (25, 50, 100, 200, 400):
    top = order[:k]
    hit = sum(1 for g in top if g in setA)
    out.append('  top-%3d out-degree: %3d/%3d = %.1f%% both-up' % (k, hit, k, 100 * hit / k))

# in-degree counterpart
ind = np.array([in_deg.get(g, 0) for g in meas], float)
rho2, p2_ = stats.spearmanr(ind, neglog)
out.append('\nSpearman(in-degree, -log10 Fisher p) = %.3f  p=%.2e' % (rho2, p2_))
order2 = sorted(meas, key=lambda g: -in_deg.get(g, 0))
for k in (50, 100, 200):
    hit = sum(1 for g in order2[:k] if g in setA)
    out.append('  top-%3d in-degree : %3d/%3d = %.1f%% both-up' % (k, hit, k, 100 * hit / k))

nine = ['STAT1', 'FOSL1', 'CXCL8', 'MMP9', 'E2F1', 'CXCL10', 'JUN', 'RELA', 'ETS1']
out.append('\nFisher top-9 out-degree:')
for g in nine:
    rk = order.index(g) + 1 if g in order else None
    out.append('  %-7s out=%4d  rank=%s/%d  in=%d' % (g, out_deg.get(g, 0), rk, len(order), in_deg.get(g, 0)))

# how many of the top-12 out-degree hubs are in the both-up set?
out.append('\ntop-12 out-degree hubs -> both-up?')
for g in order[:12]:
    out.append('  %-8s out=%4d  both-up=%s  FisherP=%.1e' % (g, out_deg.get(g, 0), g in setA, fisher[g]))

open(r'./_deg_enrich.txt', 'w', encoding='utf-8').write('\n'.join(out))
print('\n'.join(out))
