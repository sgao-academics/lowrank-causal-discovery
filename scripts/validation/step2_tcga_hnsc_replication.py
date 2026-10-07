# -*- coding: utf-8 -*-
"""External validation #1: TCGA-HNSC (n=566, geographically independent).
Tumour vs adjacent-normal for the inflammatory core, and replication of the
334 'significant in both Indian cohorts' genes."""
import os, sys, json, csv
import numpy as np
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')
T  = r'{DATA_ROOT}/cancer_application/data/TCGA_HNSC_HiSeqV2.tsv'
GW = r'{DATA_ROOT}/cancer_application/data/indian_oral/indian_genomewide_de.json'
EDGE = r'./checkpoints/edge_list.csv'
out=[]

# --- load TCGA-HNSC (all genes) ---
hdr=open(T,encoding='utf-8',errors='ignore').readline().rstrip('\n').split('\t')
samples=hdr[1:]
is_tum=np.array([s[13:15]=='01' for s in samples]); is_norm=np.array([s[13:15]=='11' for s in samples])
out.append('TCGA-HNSC samples=%d tumor=%d normal=%d'%(len(samples),is_tum.sum(),is_norm.sum()))

E={}
with open(T,encoding='utf-8',errors='ignore') as f:
    f.readline()
    for line in f:
        p=line.rstrip('\n').split('\t')
        if len(p)!=(len(samples)+1): continue
        E[p[0]]=np.array(p[1:],dtype=float)
out.append('genes loaded = %d'%len(E))

def mwu(g):
    x=E[g][is_tum]; y=E[g][is_norm]
    u,p=stats.mannwhitneyu(x,y,alternative='two-sided')
    return x.mean()-y.mean(), p

HL=['STAT1','JUN','FOSL1','RELA','CXCL8','CXCL10','MMP9','PTGS2','E2F1','ETS1','NFKB1','STAT3','SP1','AR']
out.append('\n=== inflammatory-core genes: tumour vs normal (TCGA-HNSC) ===')
out.append('%-8s %9s %11s'%('gene','delta','p'))
for g in HL:
    if g in E:
        d,p=mwu(g); out.append('%-8s %9.2f %11.2e'%(g,d,p))
    else: out.append('%-8s %9s %11s'%(g,'NA','NA'))

# --- Indian both-up set ---
gw=json.load(open(GW,encoding='utf-8'))
P1=gw['p']['GSE85195']; P2=gw['p']['GSE23558']; D1=gw['delta']['GSE85195']; D2=gw['delta']['GSE23558']
net=[]
with open(EDGE,encoding='utf-8',errors='ignore') as f:
    for row in csv.DictReader(f):
        net.append(row['source_gene']); net.append(row['target_gene'])
net=sorted(set(net))
def bh(pv,genes):
    g=np.array([pv[x] for x in genes]); o=np.argsort(g); m=len(g); q=np.empty(m); prev=1.0
    for i in range(m-1,-1,-1):
        idx=o[i]; prev=min(prev,g[idx]*m/(i+1)); q[idx]=prev
    return {genes[i]:q[i] for i in range(m)}
meas=[x for x in net if x in P1 and x in P2]
q1=bh(P1,meas); q2=bh(P2,meas)
both_up=[x for x in meas if q1[x]<0.05 and q2[x]<0.05 and D1[x]>0 and D2[x]>0]
out.append('\nIndian both-up (FDR<0.05) = %d genes'%len(both_up))

rep=[]; tested=0; up=0; sig=0
for g in both_up:
    if g not in E: continue
    tested+=1; d,p=mwu(g)
    if d>0: up+=1
    if d>0 and p<0.05: sig+=1
    rep.append((g,d,p))
out.append('measurable in TCGA-HNSC = %d'%tested)
out.append('  up (delta>0)        = %d (%.1f%%)'%(up,100*up/max(tested,1)))
out.append('  up AND p<0.05       = %d (%.1f%%)'%(sig,100*sig/max(tested,1)))

# background over same number of random genes
rng=np.random.default_rng(0)
bgs=set(rng.choice(list(E.keys()), size=min(4000,len(E)), replace=False))
bgup=0; bgsig=0; bgn=0
for g in bgs:
    v=E[g]
    x=v[is_tum]; y=v[is_norm]
    if len(x)<3 or len(y)<3: continue
    bgn+=1
    if x.mean()-y.mean()>0: bgup+=1
    u,p=stats.mannwhitneyu(x,y,alternative='two-sided')
    if x.mean()-y.mean()>0 and p<0.05: bgsig+=1
out.append('background (random %d genes): up=%.1f%%  up&p<0.05=%.1f%%'%(bgn,100*bgup/bgn,100*bgsig/bgn))

# --- 2x2 for odds ratio of both-up genes replicating ---
a=sig; b=tested-sig; c=bgsig; d=bgn-bgsig
OR=(a*d)/(b*c) if b*c>0 else float('nan')
out.append('odds ratio (both-up replicates vs background) = %.2f'%OR)

# core score
core=['STAT1','JUN','FOSL1','RELA','CXCL8','CXCL10','MMP9','PTGS2','E2F1','ETS1']
pres=[g for g in core if g in E]
z=np.vstack([(E[g]-E[g].mean())/E[g].std() for g in pres]); sc=z.mean(0)
u,p=stats.mannwhitneyu(sc[is_tum],sc[is_norm],alternative='two-sided')
out.append('\ncore-score (%d genes): tumour %.3f vs normal %.3f  p=%.2e'%(len(pres),sc[is_tum].mean(),sc[is_norm].mean(),p))

open(r'./_tcga_hnsc.txt','w',encoding='utf-8').write('\n'.join(out))
print('\n'.join(out))
