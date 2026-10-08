# -*- coding: utf-8 -*-
"""Section S5 of the supplementary material: the Benjamini-Hochberg correction
over the whole transcriptome, reported there for transparency next to the
network-scoped correction used in the main text.

Reproduces: 9,330 of 19,565 genes in GSE85195 and 4,066 of 19,565 in GSE23558
at a transcriptome-wide FDR of 5%, and the per-gene transcriptome-wide FDR of
every gene the main text emphasises.
"""
import os, sys, json
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')

GW = r'{DATA_ROOT}/cancer_application/data/indian_oral/indian_genomewide_de.json'
OUT = r'./validation_outputs/step12_transcriptome_wide_correction.txt'

FOCUS = {
    'GSE85195': ['FOSL1', 'RELA', 'CXCL8', 'MMP9', 'STAT1', 'CXCL10', 'JUN',
                 'PTGS2', 'ETS1'],
    'GSE23558': ['CA9', 'MMP1', 'SPP1', 'STAT1', 'CD274', 'MMP3', 'E2F1', 'JUN',
                 'FOSL1', 'CXCL10', 'ETS1', 'VEGFA', 'MMP10', 'CXCL8', 'PTGS2',
                 'MMP9', 'RELA'],
}


def bh(pvals):
    """Benjamini-Hochberg adjusted p-values, same convention as the paper."""
    a = np.asarray(pvals, dtype=float)
    m = a.size
    order = np.argsort(a)
    q = np.empty(m, dtype=float)
    prev = 1.0
    for i in range(m - 1, -1, -1):
        idx = order[i]
        prev = min(prev, a[idx] * m / (i + 1))
        q[idx] = prev
    return q


gw = json.load(open(GW, encoding='utf-8'))
P, D = gw['p'], gw['delta']
out = []

for acc in ('GSE85195', 'GSE23558'):
    p, d = P[acc], D[acc]
    genes = list(p.keys())
    q = dict(zip(genes, bh([p[g] for g in genes])))
    n_sig = sum(1 for g in genes if q[g] < 0.05)
    n_up = sum(1 for g in genes if q[g] < 0.05 and d[g] > 0)
    out.append('=== %s ===' % acc)
    out.append('  genes tested                    %d' % len(genes))
    out.append('  transcriptome-wide FDR < 5%%    %d  (up %d, down %d)'
               % (n_sig, n_up, n_sig - n_up))
    out.append('  per-gene transcriptome-wide FDR:')
    for g in FOCUS[acc]:
        if g in q:
            out.append('    %-7s FDR = %.3g   (delta %+.2f)' % (g, q[g], d[g]))
        else:
            out.append('    %-7s absent from this matrix' % g)
    out.append('')

txt = '\n'.join(out)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
open(OUT, 'w', encoding='utf-8').write(txt)
print(txt)
