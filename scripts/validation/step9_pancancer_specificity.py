# -*- coding: utf-8 -*-
"""External validation #4 (Section S12): is the inflammatory core oral-specific?

The 334-gene set was defined on Indian oral-cavity tumours only, so it is an
unbiased probe of other tumour types.  The same tumour-versus-adjacent-normal test
is applied to every TCGA cohort in the UCSC Xena release that has at least three
adjacent normals.  Replication = fraction of the measurable core genes that move in
the same direction and reach p < 0.05, with the odds ratio taken against an equal
number of genes drawn at random from the same cohort and platform.

Input : {DATA_ROOT}/cancer_application/data/TCGA_<CODE>_HiSeqV2.tsv (33 cohorts)
        {DATA_ROOT}/cancer_application/data/indian_oral/indian_genomewide_de.json
        ./checkpoints/edge_list.csv
Output: ./validation_outputs/step9_pancancer_specificity.txt / .json
"""
import os, sys, json, csv, math, time
import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
DATA = r'{DATA_ROOT}/cancer_application/data'
EDGE = r'./checkpoints/edge_list.csv'
GW = r'{DATA_ROOT}/cancer_application/data/indian_oral/indian_genomewide_de.json'
OUT = r'./validation_outputs/step9_pancancer_specificity.txt'
OUTJ = r'./validation_outputs/step9_pancancer_specificity.json'

ALL_33 = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC',
          'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV',
          'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM',
          'UCEC', 'UCS', 'UVM']

FISHER9 = ['STAT1', 'FOSL1', 'CXCL8', 'MMP9', 'E2F1', 'CXCL10', 'JUN', 'RELA', 'ETS1']

L = []
def say(s=''):
    L.append(s); print(s)


def load_matrix(path):
    """UCSC-Xena HiSeqV2: rows = genes, columns = samples."""
    with open(path, encoding='utf-8', errors='ignore') as f:
        hdr = f.readline().rstrip('\n').split('\t')
        samples = hdr[1:]
        genes, mats = [], []
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            i = line.find('\t')
            if i <= 0:
                continue
            v = np.fromstring(line[i + 1:], sep='\t', dtype=np.float32)
            if v.size != len(samples):
                continue
            genes.append(line[:i]); mats.append(v)
    return genes, (np.vstack(mats) if mats else np.zeros((0, len(samples)), np.float32)), samples


net = set()
with open(EDGE, encoding='utf-8', errors='ignore') as f:
    for row in csv.DictReader(f):
        net.add(row['source_gene']); net.add(row['target_gene'])

gw = json.load(open(GW, encoding='utf-8'))
P1, P2 = gw['p']['GSE85195'], gw['p']['GSE23558']
D1, D2 = gw['delta']['GSE85195'], gw['delta']['GSE23558']


def bh(pv, genes):
    g = np.array([pv[x] for x in genes]); o = np.argsort(g); m = len(g)
    q = np.empty(m); prev = 1.0
    for i in range(m - 1, -1, -1):
        idx = o[i]
        prev = min(prev, g[idx] * m / (i + 1))
        q[idx] = prev
    return {genes[i]: q[i] for i in range(m)}


meas = [x for x in sorted(net) if x in P1 and x in P2]
q1, q2 = bh(P1, meas), bh(P2, meas)
both_up = [x for x in meas if q1[x] < 0.05 and q2[x] < 0.05 and D1[x] > 0 and D2[x] > 0]
both_sig = [x for x in meas if q1[x] < 0.05 and q2[x] < 0.05]
say('network genes = %d ; measurable in Indian = %d' % (len(net), len(meas)))
say('Indian both-up FDR<0.05 (334-set) = %d ; both-significant (578-set) = %d'
    % (len(both_up), len(both_sig)))
say('Fisher Top-9 = %s' % ', '.join(FISHER9))
say('')

CORE = set(both_up) | set(FISHER9)
rng = np.random.default_rng(20261007)

res = {}
t_all = time.time()
for code in ALL_33:
    p = os.path.join(DATA, 'TCGA_%s_HiSeqV2.tsv' % code)
    if not os.path.exists(p):
        say('%-5s  MISSING FILE' % code); continue
    t0 = time.time()
    genes, M, samples = load_matrix(p)
    idx = {g: i for i, g in enumerate(genes)}
    is_t = np.array([s[13:15] == '01' for s in samples])
    is_n = np.array([s[13:15] == '11' for s in samples])
    nt, nn = int(is_t.sum()), int(is_n.sum())
    if nt < 5 or nn < 3:
        say('%-5s  tumour=%d normal=%d  -> SKIPPED (too few normals)' % (code, nt, nn))
        res[code] = {'tumor': nt, 'normal': nn, 'skipped': True}
        continue

    Xt, Xn = M[:, is_t], M[:, is_n]

    def scan(gi_list):
        up = sig = 0; per = {}
        for g in gi_list:
            x = Xt[idx[g]]; y = Xn[idx[g]]
            d = float(x.mean() - y.mean())
            try:
                _, pv = stats.mannwhitneyu(x, y, alternative='two-sided')
            except ValueError:
                pv = 1.0
            pv = float(pv)
            per[g] = [round(d, 4), pv]
            if d > 0:
                up += 1
            if d > 0 and pv < 0.05:
                sig += 1
        return up, sig, per

    core_meas = [g for g in sorted(CORE) if g in idx]
    up, sig, per = scan(core_meas)

    allg = [g for g in genes if g not in CORE]
    bg = list(rng.choice(allg, size=min(len(core_meas), len(allg)), replace=False))
    bup, bsig, _ = scan(bg)
    btest = len(bg)

    n_core = len(core_meas)
    a_, b_ = sig, n_core - sig
    c_, d_ = bsig, btest - bsig
    OR = (a_ * d_) / (b_ * c_) if b_ * c_ > 0 else float('nan')

    pres = [g for g in FISHER9 if g in idx]
    if len(pres) >= 5:
        Z = np.vstack([(M[idx[g]] - M[idx[g]].mean()) / (M[idx[g]].std() + 1e-9) for g in pres])
        sc = Z.mean(0)
        _, pcs = stats.mannwhitneyu(sc[is_t], sc[is_n], alternative='two-sided')
        cs_d = float(sc[is_t].mean() - sc[is_n].mean()); pcs = float(pcs)
    else:
        pcs, cs_d = float('nan'), float('nan')

    res[code] = dict(tumor=nt, normal=nn, n_core=n_core, n_up=up, n_sig=sig,
                     bg_test=btest, bg_up=bup, bg_sig=bsig,
                     OR=None if math.isnan(OR) else round(OR, 2),
                     core_score_delta=None if math.isnan(cs_d) else round(cs_d, 3),
                     core_score_p=None if math.isnan(pcs) else pcs, per_gene=per)
    say('%-5s n=%3d/%3d | core %3d  up %3d(%3.0f%%)  up&p<.05 %3d(%3.0f%%) | bg %.0f%%/%.0f%% (n=%d) | OR %5s | coreScore d=%+6.3f p=%.2e | %.0fs'
        % (code, nt, nn, n_core, up, 100 * up / max(n_core, 1), sig,
           100 * sig / max(n_core, 1), 100 * bup / max(btest, 1), 100 * bsig / max(btest, 1),
           btest, ('%.2f' % OR) if not math.isnan(OR) else 'NA',
           0 if math.isnan(cs_d) else cs_d, 0 if math.isnan(pcs) else pcs, time.time() - t0))

say('')
say('total elapsed %.1f s' % (time.time() - t_all))
say('')
say('=== ranking by replication strength (fraction of core genes up&p<0.05) ===')
rank = [(c, r['n_sig'] / max(r['n_core'], 1), r['OR'], r['tumor'], r['normal'])
        for c, r in res.items() if not r.get('skipped') and r.get('n_core')]
rank.sort(key=lambda x: -x[1])
for i, (c, fr, OR, nt, nn) in enumerate(rank, 1):
    say('%2d. %-5s  %3.0f%%  OR=%s   (n=%d/%d)' % (i, c, 100 * fr, OR, nt, nn))
say('')
say('=== ranking by core-score delta (Fisher-9, tumour vs normal) ===')
rank2 = [(c, r['core_score_delta'], r['core_score_p'])
         for c, r in res.items() if not r.get('skipped') and r.get('core_score_delta') is not None]
rank2.sort(key=lambda x: -x[1])
for i, (c, d, p) in enumerate(rank2, 1):
    say('%2d. %-5s  d=%+.3f  p=%.2e' % (i, c, d, p))

os.makedirs('./validation_outputs', exist_ok=True)
json.dump({'both_up_334': both_up, 'both_sig_578': both_sig, 'fisher9': FISHER9,
           'cohorts': res}, open(OUTJ, 'w', encoding='utf-8'), indent=1)
open(OUT, 'w', encoding='utf-8').write('\n'.join(L))
print('\nwritten', OUT)
