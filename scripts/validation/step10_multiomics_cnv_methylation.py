# -*- coding: utf-8 -*-
"""External validation #5 (Section S13): copy number and DNA methylation of the core.

The core is defined on mRNA alone, and NFKB1 mRNA does not move while its targets
do.  A reviewer will ask whether the activation could instead be a copy-number
gain or a promoter demethylation of the regulators.  Both layers are tested in
TCGA-HNSC.

(A) CNV      -- GISTIC thresholded calls per gene; frequency of amplification and
                deletion for the core genes against the genome-wide background of
                the same file.
(B) methylation -- promoter probes (TSS1500, TSS200, first exon, 5'UTR) taken from
                the Illumina HumanMethylation450 manifest, tested tumour versus
                adjacent normal (two-sided Mann-Whitney).

Both source files are downloaded on first use from the UCSC Xena HNSC sample map
and from Illumina; nothing is redistributed.
"""
import os, sys, gzip, csv, json, time, urllib.request
import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
DATA = r'{DATA_ROOT}/cancer_application/data/hnsc_multiomics'
OUT = r'./validation_outputs/step10_multiomics_cnv_methylation.txt'
OUTJ = r'./validation_outputs/step10_multiomics_cnv_methylation.json'
XENA = 'https://tcga.xenahubs.net/download/TCGA.HNSC.sampleMap/'
MAN_URL = ('https://webdata.illumina.com/downloads/productfiles/humanmethylation450/'
           'humanmethylation450_15017482_v1-2.csv')

CNV_GZ = 'HNSC_CNV_gistic.by_genes.gz'
MET_GZ = 'HNSC_methylation450.gz'
MAN = 'hm450_manifest.csv'

CORE9 = ['STAT1', 'FOSL1', 'CXCL8', 'MMP9', 'E2F1', 'CXCL10', 'JUN', 'RELA', 'ETS1']
REG = ['NFKB1', 'NFKB2', 'RELA', 'STAT1', 'STAT3', 'JUN', 'FOSL1', 'ETS1', 'E2F1',
       'SP1', 'MYC', 'TP53']
PRO = ('TSS1500', 'TSS200', '1stExon', "5'UTR")
# methylation arm: the nine core genes plus the NF-kB regulators they act through
METH_GENES = CORE9 + ['NFKB1', 'NFKB2', 'STAT3']

L = []
def say(s=''):
    L.append(s); print(s)


def fetch(url, path, min_size=1000):
    if os.path.exists(path) and os.path.getsize(path) > min_size:
        return path
    os.makedirs(os.path.dirname(path), exist_ok=True)
    print('downloading', os.path.basename(path), '...')
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    try:
        with urllib.request.urlopen(req, timeout=300) as r, open(path, 'wb') as f:
            n = 0
            while True:
                b = r.read(1 << 20)
                if not b:
                    break
                f.write(b); n += len(b)
    except Exception as e:
        raise SystemExit('%s could not be downloaded (%s).\nRetrieve it manually from %s\n'
                         'and place it at %s, then re-run.' % (os.path.basename(path), e, url, path))
    print('  %.1f MB' % (n / 1048576))
    return path


t0 = time.time()
os.makedirs(DATA, exist_ok=True)
fetch(XENA + CNV_GZ, os.path.join(DATA, CNV_GZ))

# ---------------------------------------------------------------- (A) CNV
say('=== (A) GISTIC copy-number, TCGA-HNSC ===')
with gzip.open(os.path.join(DATA, CNV_GZ), 'rt', encoding='utf-8', errors='ignore') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    samples = hdr[1:]
    tgt = set(REG) | set(CORE9)
    rows = {}
    amp = dele = ncall = 0
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) != len(samples) + 1:
            continue
        v = np.fromstring('\t'.join(p[1:]), sep='\t', dtype=np.float32)
        if v.size != len(samples):
            continue
        a = int((v >= 1).sum()); d = int((v <= -1).sum())
        amp += a; dele += d; ncall += v.size
        if p[0] in tgt:
            rows[p[0]] = (a, d, v.size)
say('samples = %d ; genes in file = %d' % (len(samples), ncall // max(len(samples), 1)))
say('genome-wide: amplified calls %.3f%%  deleted calls %.3f%%' % (100 * amp / ncall, 100 * dele / ncall))
say('%-8s %8s %8s %8s' % ('gene', 'amp%', 'del%', 'n'))
cnv_out = {}
for g in sorted(set(REG) | set(CORE9)):
    if g in rows:
        a, d, n = rows[g]
        cnv_out[g] = {'amp_pct': round(100 * a / n, 2), 'del_pct': round(100 * d / n, 2), 'n': n}
        say('%-8s %8.2f %8.2f %8d' % (g, 100 * a / n, 100 * d / n, n))
    else:
        say('%-8s %8s %8s' % (g, 'NA', 'NA'))
say('')

# ---------------------------------------------------------------- (B) methylation
say('=== (B) 450k promoter methylation of the core genes, TCGA-HNSC ===')
fetch(MAN_URL, os.path.join(DATA, MAN))
probes = {}
with open(os.path.join(DATA, MAN), encoding='utf-8', errors='ignore', newline='') as f:
    rd = csv.reader(f)
    hdr = None
    for row in rd:
        if row and row[0] == 'IlmnID':
            hdr = row; break
    ia = hdr.index('IlmnID'); ig = hdr.index('UCSC_RefGene_Name'); igp = hdr.index('UCSC_RefGene_Group')
    for row in rd:
        if len(row) <= max(ia, ig, igp):
            continue
        for g, grp in zip(row[ig].split(';'), row[igp].split(';')):
            if g in METH_GENES and grp in PRO:
                probes.setdefault(g, set()).add(row[ia])
say('promoter probes: %s' % ', '.join('%s=%d' % (g, len(probes.get(g, ()))) for g in METH_GENES))
allp = set()
for v in probes.values():
    allp |= v

fetch(XENA + MET_GZ, os.path.join(DATA, MET_GZ))
me = {}
with gzip.open(os.path.join(DATA, MET_GZ), 'rt', encoding='utf-8', errors='ignore') as f:
    h = f.readline().rstrip('\n').split('\t')
    msamples = h[1:]
    for line in f:
        i = line.find('\t')
        if i <= 0:
            continue
        pid = line[:i]
        if pid in allp:
            raw = line[i + 1:].rstrip('\n').split('\t')
            if len(raw) != len(msamples):
                continue
            try:
                me[pid] = np.array([np.nan if t in ('', 'NA', 'NaN') else float(t) for t in raw],
                                   np.float32)
            except ValueError:
                continue
is_t = np.array([s[13:15] == '01' for s in msamples])
is_n = np.array([s[13:15] == '11' for s in msamples])
say('methylation samples %d (tumour %d / normal %d) ; probes found %d / %d'
    % (len(msamples), is_t.sum(), is_n.sum(), len(me), len(allp)))

me_out = {}
say('')
say('%-7s %-12s %8s %8s %9s %10s' % ('gene', 'probe', 'beta_T', 'beta_N', 'delta', 'p'))
for g in METH_GENES:
    recs = []
    for p in sorted(probes.get(g, ())):
        if p not in me:
            continue
        v = me[p]; xt = v[is_t]; xn = v[is_n]
        xt = xt[~np.isnan(xt)]; xn = xn[~np.isnan(xn)]
        if xt.size < 3 or xn.size < 3:
            continue
        bt, bn = float(xt.mean()), float(xn.mean())
        _, pv = stats.mannwhitneyu(xt, xn, alternative='two-sided')
        recs.append({'probe': p, 'beta_tumour': round(bt, 4), 'beta_normal': round(bn, 4),
                     'delta': round(bt - bn, 4), 'p': float(pv)})
        say('%-7s %-12s %8.3f %8.3f %+9.3f %10.2e' % (g, p, bt, bn, bt - bn, pv))
    if recs:
        d = [x['delta'] for x in recs]
        me_out[g] = {'n_probes': len(recs), 'delta_min': min(d), 'delta_max': max(d),
                     'delta_mean': round(float(np.mean(d)), 4),
                     'n_nominal_p05': int(sum(1 for x in recs if x['p'] < 0.05)), 'probes': recs}
        say('   -> %-6s n=%2d  delta range %+.4f .. %+.4f  mean %+.4f'
            % (g, len(recs), min(d), max(d), float(np.mean(d))))

say('')
say('elapsed %.1f s' % (time.time() - t0))
os.makedirs('./validation_outputs', exist_ok=True)
json.dump({'cnv': cnv_out, 'methylation': me_out,
           'genome_wide_amp_pct': round(100 * amp / ncall, 3),
           'genome_wide_del_pct': round(100 * dele / ncall, 3),
           'methylation_samples': {'tumour': int(is_t.sum()), 'normal': int(is_n.sum())}},
          open(OUTJ, 'w', encoding='utf-8'), indent=1)
open(OUT, 'w', encoding='utf-8').write('\n'.join(L))
print('written', OUT)
