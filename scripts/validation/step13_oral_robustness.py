# -*- coding: utf-8 -*-
"""Section S4 of the supplementary material: robustness of the Indian-cohort
analysis.

Reproduces, from the primary data:

  * the sex composition of the two cohorts and of the groups being compared,
    read from the GEO series-matrix characteristics;
  * the Y-chromosome check, as a single score over all chrY probes rather than
    gene by gene, since the individual transcripts are not individually
    informative at n = 5 controls;
  * the expression-level matching of the network against the remaining
    transcriptome, within deciles of mean expression;
  * the out-degree enrichment for up-regulation after every X-, Y- and
    mitochondrial gene has been removed.

The up-regulation criterion is the one used throughout Section S4: a raw
two-sided Mann-Whitney p < 0.05 together with a positive difference of group
means. Chromosome assignment comes from the UCSC hg19 refGene table; the
platform's cyto-band column cannot be parsed into a chromosome reliably.

The null for the enrichment is a binomial draw of k genes from the
sex-chromosome-excluded network pool, at that pool's own up-regulation rate, so
a single background rate is used for every cutoff.
"""
import os, sys, gzip, csv, json, collections
import numpy as np
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')

DATA = r'{DATA_ROOT}/cancer_application/data/indian_oral'
REFGENE = r'{DATA_ROOT}/cancer_application/data/hg19_refGene.txt.gz'
ANNOT = os.path.join(DATA, 'GPL6480.annot.gz')
GW = os.path.join(DATA, 'indian_genomewide_de.json')
EDGE = r'./checkpoints/edge_list.csv'
OUTDIR = r'./validation_outputs'
SERIES = {'GSE85195': os.path.join(DATA, 'GSE85195_series_matrix.txt.gz'),
          'GSE23558': os.path.join(DATA, 'GSE23558_series_matrix.txt.gz')}

SEXCHROMS = {'chrX', 'chrY', 'chrM', 'chrMT'}
Y_MARKERS = ['RPS4Y1', 'DDX3Y', 'KDM5D', 'UTY', 'USP9Y', 'EIF1AY', 'ZFY', 'NLGN4Y']


def gpl_symbols(path):
    """probe name -> gene symbol, from the GEO platform annotation.

    The symbol is taken verbatim, as deposited: the differential-expression
    tables this section is built on were computed on this same convention, so
    splitting multi-symbol entries would compare against a different gene
    universe.
    """
    sym, in_tab, idx_sym = {}, False, None
    with gzip.open(path, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('!platform_table_begin'):
                in_tab = True
                continue
            if line.startswith('!platform_table_end'):
                break
            if not in_tab:
                continue
            c = line.split('\t')
            if idx_sym is None:
                low = [x.strip().lower() for x in c]
                for cand in ('gene_symbol', 'gene symbol', 'symbol'):
                    if cand in low:
                        idx_sym = low.index(cand)
                        break
                continue
            if idx_sym is None or len(c) <= idx_sym:
                continue
            g = c[idx_sym].strip().upper()
            if g and g not in ('', '---'):
                sym.setdefault(c[0].strip(), g)
    return sym


def chrom_map(path):
    """gene symbol -> chromosome, from the UCSC refGene table."""
    m = {}
    with gzip.open(path, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            c = line.rstrip('\n').split('\t')
            if len(c) >= 13 and c[12] and c[12] not in m:
                m[c[12]] = c[2]
    return m


def series_matrix(path):
    """sample titles, source names, sex labels and the probe x sample matrix."""
    titles, src, sex, vals, in_tab = [], [], [], {}, False
    with gzip.open(path, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('!Sample_title'):
                titles = [x.strip().strip('"') for x in line.split('\t')[1:]]
            elif line.startswith('!Sample_source_name'):
                src = [x.strip().strip('"') for x in line.split('\t')[1:]]
            elif line.startswith('!Sample_characteristics'):
                row = [x.strip().strip('"') for x in line.split('\t')[1:]]
                if row and row[0].lower().startswith(('gender', 'sex')):
                    sex = [x.split(':', 1)[-1].strip() for x in row]
            if line.startswith('!series_matrix_table_begin'):
                in_tab = True
                continue
            if line.startswith('!series_matrix_table_end'):
                break
            if not in_tab:
                continue
            c = line.split('\t')
            if 'header' not in vals:
                vals['header'] = True
                continue
            pid = c[0].strip().strip('"')
            try:
                vals[pid] = [float(x) for x in c[1:1 + len(titles)]]
            except ValueError:
                continue
    vals.pop('header', None)
    return titles, src, sex, vals


def gene_means(probes, sym):
    """expression per gene symbol; where several probes map to one symbol the
    probe with the highest mean signal is retained, as in Methods."""
    best = {}
    for pid, v in probes.items():
        g = sym.get(pid)
        if not g:
            continue
        m = float(np.mean(v))
        if g not in best or m > best[g][0]:
            best[g] = (m, v)
    return {g: b[1] for g, b in best.items()}


sym = gpl_symbols(ANNOT)
cm = chrom_map(REFGENE)
gw = json.load(open(GW, encoding='utf-8'))
P, D = gw['p'], gw['delta']

net = set()
out_edge = set()          # genes with at least one outgoing edge of any type
od = collections.Counter()
for r in csv.DictReader(open(EDGE, encoding='utf-8', errors='ignore')):
    net.add(r['source_gene'])
    net.add(r['target_gene'])
    out_edge.add(r['source_gene'])
    if r['source'] == 'discovered':
        od[r['source_gene']] += 1

out = []
report = {}
out.append('gene symbols on the platform: %d' % len(sym))
out.append('network genes in the released edge list: %d' % len(net))

# ------------------------------------------------------------ sex composition
out.append('')
out.append('=== sex composition (GEO series-matrix characteristics) ===')
comp = {}
for acc in ('GSE85195', 'GSE23558'):
    titles, src, sex, _ = series_matrix(SERIES[acc])
    groups = collections.OrderedDict()
    for s in dict.fromkeys(src):
        groups[s] = collections.Counter(
            sex[i] if i < len(sex) else 'n/a'
            for i in range(len(src)) if src[i] == s)
    comp[acc] = {g: dict(c) for g, c in groups.items()}
    out.append('%s (n = %d)' % (acc, len(src)))
    for g, c in groups.items():
        out.append('   %-40s %s' % (g[:40], dict(c)))
report['sex_composition'] = comp

# ------------------------------------------------------------------ Y check
out.append('')
out.append('=== Y-chromosome check, GSE23558 ===')
titles, src, sex, probes = series_matrix(SERIES['GSE23558'])
tum = [i for i in range(len(src)) if 'umor' in src[i]]
ctl = [i for i in range(len(src)) if 'ormal' in src[i]]
yprobes = [pid for pid in probes if cm.get(sym.get(pid, ''), '') == 'chrY']
out.append('probes mapping to a chrY gene: %d over %d distinct genes'
           % (len(yprobes), len({sym[pid] for pid in yprobes if sym.get(pid)})))
for g in Y_MARKERS:
    if g in P['GSE23558']:
        out.append('   %-8s delta %+6.2f  p = %.3g'
                   % (g, D['GSE23558'][g], P['GSE23558'][g]))
score = np.array([probes[pid] for pid in yprobes], dtype=float).mean(axis=0)
u, pv = stats.mannwhitneyu(score[tum], score[ctl], alternative='two-sided')
out.append('combined Y score (mean of %d probes): tumour median %.3f, '
           'control median %.3f, Mann-Whitney p = %.3g'
           % (len(yprobes), float(np.median(score[tum])),
              float(np.median(score[ctl])), pv))
report['y_check'] = {'n_probes': len(yprobes),
                     'tumour_median': float(np.median(score[tum])),
                     'control_median': float(np.median(score[ctl])),
                     'p': float(pv)}

# -------------------------------------------------------- expression deciles
out.append('')
out.append('=== expression-level matching, within deciles of mean expression ===')
deciles = {}
for acc in ('GSE85195', 'GSE23558'):
    _, _, _, probes = series_matrix(SERIES[acc])
    gm = gene_means(probes, sym)
    pv, dl = P[acc], D[acc]
    meas = [g for g in gm if g in pv and g in dl]
    # the decile comparison follows the convention of the deposited table: the
    # "network" side is the set of genes that carry an outgoing edge, and the
    # comparator is every remaining measurable gene (which includes the genes
    # that appear in the network only as targets, making the test conservative).
    netm = [g for g in out_edge if g in meas]
    oth = [g for g in meas if g not in out_edge]
    means = {g: float(np.mean(gm[g])) for g in meas}
    allg = netm + oth
    qs = np.quantile([means[g] for g in allg], np.linspace(0, 1, 11))
    rows = []
    out.append('%s (%d measurable genes)' % (acc, len(meas)))
    out.append('   %-8s %7s %7s %7s %7s %7s'
               % ('decile', 'n_net', 'n_oth', 'up_net', 'up_oth', 'ratio'))
    for k in range(10):
        lo, hi = qs[k], qs[k + 1]
        a = [g for g in netm if lo <= means[g] <= hi]
        b = [g for g in oth if lo <= means[g] <= hi]
        if len(a) < 10 or len(b) < 10:
            continue
        ua = float(np.mean([1.0 if (pv[g] < 0.05 and dl[g] > 0) else 0.0 for g in a]))
        ub = float(np.mean([1.0 if (pv[g] < 0.05 and dl[g] > 0) else 0.0 for g in b]))
        rr = round(ua / ub, 2) if ub > 0 else None
        rows.append({'decile': k + 1, 'n_net': len(a), 'n_oth': len(b),
                     'up_net': round(ua, 3), 'up_oth': round(ub, 3), 'ratio': rr})
        out.append('   %-8d %7d %7d %7.3f %7.3f %7s'
                   % (k + 1, len(a), len(b), ua, ub, rr))
    deciles[acc] = rows
report['expr_matched'] = deciles

# --------------------------------------------- out-degree enrichment, no sex
out.append('')
out.append('=== out-degree enrichment for up-regulation, X/Y/MT removed ===')
enr = {}
for acc in ('GSE85195', 'GSE23558'):
    pv, dl = P[acc], D[acc]
    meas = [g for g in net if g in pv and g in dl]
    keep = [g for g in meas if cm.get(g, '') not in SEXCHROMS]
    n_sex = len(meas) - len(keep)
    up = lambda g: 1.0 if (dl[g] > 0 and pv[g] < 0.05) else 0.0
    bg = float(np.mean([up(g) for g in keep]))       # one background for all k
    order = sorted([g for g in keep if od.get(g, 0) > 0], key=lambda g: -od[g])
    res = {}
    out.append('%s: measurable network genes %d -> %d (%d X/Y/MT removed)'
               % (acc, len(meas), len(keep), n_sex))
    out.append('   network background (X/Y/MT removed) = %.4f' % bg)
    for k in (50, 100, 200):
        top = order[:k]
        obs = float(np.mean([up(g) for g in top]))
        hits = int(round(obs * k))
        pb = stats.binomtest(hits, k, bg, alternative='greater').pvalue
        res[str(k)] = {'obs': round(obs, 4), 'bg': round(bg, 4),
                       'hits': hits, 'p': float(pb)}
        out.append('   top-%3d: %.3f up-regulated against %.4f (p = %.4f)'
                   % (k, obs, bg, pb))
    enr[acc] = {'n_measurable': len(meas), 'n_sex_removed': n_sex,
                'background': round(bg, 4), 'topk': res}
report['enrichment_no_sex'] = enr

os.makedirs(OUTDIR, exist_ok=True)
json.dump(report, open(os.path.join(OUTDIR, 'step13_oral_robustness.json'), 'w'),
          indent=1)
txt = '\n'.join(out)
open(os.path.join(OUTDIR, 'step13_oral_robustness.txt'), 'w',
     encoding='utf-8').write(txt)
print(txt)
