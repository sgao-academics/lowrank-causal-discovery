# -*- coding: utf-8 -*-
"""External validation #2: DepMap CRISPR (Chronos) selective dependency of the core
in head-and-neck cell lines; plus inspect the drug-repurposing matrix."""
import os, sys, csv, json
import numpy as np
import pandas as pd
sys.stdout.reconfigure(encoding='utf-8')
D = r'{DATA_ROOT}/cancer_application/data/depmap'
out=[]

# --- cell-line metadata ---
model=pd.read_csv(os.path.join(D,'Model.csv'), low_memory=False)
out.append('Model.csv shape=%s cols=%s'%(model.shape, list(model.columns)[:14]))
lin_col='OncotreeLineage'
out.append('lineages=%s'%sorted(model[lin_col].dropna().unique().tolist()))
hn=model[model[lin_col].astype(str).str.contains('Head and Neck', case=False, na=False)]
out.append('head&neck lines = %d'%len(hn))

# --- CRISPR gene effect: read header, keep needed columns ---
hdr=open(os.path.join(D,'CRISPRGeneEffect.csv'),encoding='utf-8',errors='ignore').readline().rstrip('\n').split(',')
out.append('CRISPRGeneEffect columns=%d  e.g. %s'%(len(hdr), hdr[1:4]))
# map gene symbol -> column
def sym(c): return c.split(' ')[0].strip().strip('"')
base={}
for i,c in enumerate(hdr[1:],start=1):
    base.setdefault(sym(c), i)
CORE=['STAT1','JUN','FOSL1','RELA','CXCL8','CXCL10','MMP9','PTGS2','E2F1','ETS1','NFKB1','STAT3','SP1']
cols=[0]+[base[g] for g in CORE if g in base]
use=[hdr[i] for i in cols]
out.append('core genes present in DepMap CRISPR = %d/%d'%(len(cols)-1,len(CORE)))
ge=pd.read_csv(os.path.join(D,'CRISPRGeneEffect.csv'), usecols=cols, low_memory=False)
ge.columns=['ModelID']+[sym(c) for c in use[1:]]
ge=ge.set_index('ModelID')
out.append('CRISPRGeneEffect loaded shape=%s'%str(ge.shape))

hn_ids=set(hn['ModelID'])
is_hn=ge.index.isin(hn_ids)
out.append('\n=== mean gene effect (Chronos; more negative = more dependent) ===')
out.append('%-8s %12s %12s %7s'%('gene','HN_lines','other','HN_n'))
for g in ge.columns:
    a=ge.loc[is_hn,g].dropna(); b=ge.loc[~is_hn,g].dropna()
    if len(a)<5: continue
    out.append('%-8s %12.3f %12.3f %7d'%(g,a.mean(),b.mean(),len(a)))

# how many of the 334 Indian both-up genes are selective dependencies in HN
GW=r'{DATA_ROOT}/cancer_application/data/indian_oral/indian_genomewide_de.json'
EDGE=r'./checkpoints/edge_list.csv'
gw=json.load(open(GW,encoding='utf-8'))
P1=gw['p']['GSE85195'];P2=gw['p']['GSE23558'];D1=gw['delta']['GSE85195'];D2=gw['delta']['GSE23558']
net=[]
with open(EDGE,encoding='utf-8',errors='ignore') as f:
    for row in csv.DictReader(f):
        net.append(row['source_gene']);net.append(row['target_gene'])
net=sorted(set(net))
def bh(pv,genes):
    g=np.array([pv[x] for x in genes]);o=np.argsort(g);m=len(g);q=np.empty(m);prev=1.0
    for i in range(m-1,-1,-1):
        idx=o[i];prev=min(prev,g[idx]*m/(i+1));q[idx]=prev
    return {genes[i]:q[i] for i in range(m)}
meas=[x for x in net if x in P1 and x in P2]
q1=bh(P1,meas);q2=bh(P2,meas)
both_up=[x for x in meas if q1[x]<0.05 and q2[x]<0.05 and D1[x]>0 and D2[x]>0]
out.append('\n=== are the 334 Indian both-up genes selective dependencies in H&N lines? ===')
ge_i=pd.read_csv(os.path.join(D,'CRISPRGeneEffect.csv'), low_memory=False)
ge_i.columns=['ModelID']+[sym(c) for c in hdr[1:]]
ge_i=ge_i.set_index('ModelID')
is_hn2=ge_i.index.isin(hn_ids)
present=[g for g in both_up if g in ge_i.columns]
core_m=ge_i.loc[is_hn2, present].mean(axis=0)          # HN mean effect
other_m=ge_i.loc[~is_hn2, present].mean(axis=0)
delta=(other_m-core_m)   # positive = more dependent in HN
out.append('both-up genes present in DepMap = %d'%len(present))
out.append('  mean HN gene-effect  = %.3f'%core_m.mean())
out.append('  mean other           = %.3f'%other_m.mean())
out.append('  top HN-selective (other-minus-HN most positive):')
for g in delta.sort_values(ascending=False).head(12).index:
    out.append('     %-8s HN=%.3f other=%.3f  sel=%.3f'%(g,core_m[g],other_m[g],delta[g]))
# common-essential baseline
common=pd.read_csv(os.path.join(D,'AchillesCommonEssentialControls.csv'))
out.append('  n common-essential controls = %d'%len(common))

# --- inspect drug repurposing matrix ---
rp=os.path.join(D,'REPURPOSINGLog2IC50Matrix.csv')
with open(rp,encoding='utf-8',errors='ignore') as f:
    h0=f.readline().rstrip('\n')[:300]; r1=f.readline().rstrip('\n')[:300]
out.append('\n=== REPURPOSINGLog2IC50Matrix.csv ===')
out.append('  header[:300]=%s'%h0)
out.append('  row1[:300]=%s'%r1)

open(r'./_depmap.txt','w',encoding='utf-8').write('\n'.join(out))
print('\n'.join(out)[:4000])
