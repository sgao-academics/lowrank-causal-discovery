# -*- coding: utf-8 -*-
"""Independent recomputation of the 'significant in both cohorts' gene set.
Ground truth = indian_genomewide_de.json (per-gene MWU p and delta, 19565 genes)
+ network gene set from edge_list.csv.  BH FDR<0.05 computed over the network genes."""
import json, csv, os, sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
D  = r'{DATA_ROOT}/cancer_application/data/indian_oral'
EDGE = r'./checkpoints/edge_list.csv'

out=[]
# --- 1. network genes from edge list ---
net=set(); nrow=0; srcs=[]
with open(EDGE, encoding='utf-8', errors='ignore') as f:
    r=csv.DictReader(f)
    cols=r.fieldnames
    out.append('edge_list columns = %s'%cols)
    for row in r:
        nrow+=1
        a=row.get('source_gene') or row.get('source') or row.get('from')
        b=row.get('target_gene') or row.get('target') or row.get('to')
        if a: net.add(a)
        if b: net.add(b)
out.append('edge rows=%d  unique network genes=%d'%(nrow,len(net)))

# --- 2. per-gene p/delta ---
gw=json.load(open(os.path.join(D,'indian_genomewide_de.json'),encoding='utf-8'))
P1=gw['p']['GSE85195']; P2=gw['p']['GSE23558']
D1=gw['delta']['GSE85195']; D2=gw['delta']['GSE23558']

measured=[g for g in net if g in P1 and g in P2]
out.append('network genes measurable in both = %d'%len(measured))

def bh(pvals):
    """return dict gene->q and the reject set at q<0.05"""
    g=np.array([pvals[g] for g in measured], dtype=float)
    order=np.argsort(g); m=len(g)
    q=np.empty(m); prev=1.0
    ranked=np.empty(m)
    for i in range(m-1,-1,-1):
        idx=order[i]
        val=g[idx]*m/(i+1)
        prev=min(prev,val)
        ranked[idx]=prev
    qd={measured[i]:ranked[i] for i in range(m)}
    return qd

q1=bh(P1); q2=bh(P2)
rej1={g for g in measured if q1[g]<0.05}
rej2={g for g in measured if q2[g]<0.05}
out.append('\nGSE85195: BH<0.05 = %d  (up=%d down=%d)'%(len(rej1), sum(1 for g in rej1 if D1[g]>0), sum(1 for g in rej1 if D1[g]<=0)))
out.append('GSE23558: BH<0.05 = %d  (up=%d down=%d)'%(len(rej2), sum(1 for g in rej2 if D2[g]>0), sum(1 for g in rej2 if D2[g]<=0)))

both=rej1 & rej2
both_up=sorted([g for g in both if D1[g]>0 and D2[g]>0], key=lambda g: q1[g]+q2[g])
both_dn=sorted([g for g in both if D1[g]<0 and D2[g]<0], key=lambda g: q1[g]+q2[g])
out.append('\n*** significant in BOTH cohorts: %d total ***'%len(both))
out.append('UP in both  (%d): %s'%(len(both_up), both_up))
out.append('DOWN in both (%d): %s'%(len(both_dn), both_dn))
# mixed-direction genes significant in both
mixed=[g for g in both if (D1[g]>0)!=(D2[g]>0)]
out.append('mixed/zero direction in both (%d): %s'%(len(mixed), mixed))

# --- 3. the manuscript's two stated lists ---
ABS = ['STAT1','JUN','FOSL1','RELA','CXCL8','CXCL10','MMP9','PTGS2','E2F1']   # from abstract
RES = ['STAT1','FOSL1','CXCL8','MMP9','E2F1','CXCL10','JUN','RELA','ETS1']    # from Results
out.append('\n=== check each candidate against recomputed truth ===')
out.append('%-8s %10s %10s %6s %6s  %s'%('gene','q1','q2','in1','in2','verdict'))
for tag,lst in [('ABSTRACT',ABS),('RESULTS',RES)]:
    out.append('-- %s list --'%tag)
    for g in lst:
        if g in q1 and g in q2:
            in1 = q1[g]<0.05; in2=q2[g]<0.05
            out.append('%-8s %10.2e %10.2e %6s %6s  delta=%.2f/%.2f  %s'%(
                g,q1[g],q2[g],in1,in2,D1[g],D2[g], 'BOTH' if (in1 and in2) else ('only1' if in1 else ('only2' if in2 else 'NEITHER'))))
        else:
            out.append('%-8s  NOT MEASURABLE-in-network'%g)

# --- 4. does the recomputed both_up set match either list? ---
out.append('\nABS set == recomputed both_up? %s'%(sorted(ABS)==both_up))
out.append('RES set == recomputed both_up? %s'%(sorted(RES)==both_up))

# --- 5. report q for ALL paper-highlighted genes ---
HL=['NFKB1','RELA','STAT1','STAT3','JUN','FOSL1','SP1','ETS1','E2F1','CXCL8','CXCL10','MMP9','PTGS2','IL6','CCL5','IL1B','CCND1','CDKN1A','VEGFA']
out.append('\n=== all highlighted genes, recomputed q ===')
out.append('%-8s %10s %10s %12s %12s'%('gene','q1','q2','d1','d2'))
for g in HL:
    if g in q1 and g in q2:
        out.append('%-8s %10.2e %10.2e %12.2f %12.2f'%(g,q1[g],q2[g],D1[g],D2[g]))

txt='\n'.join(out)
open(r'./_core9.txt','w',encoding='utf-8').write(txt)
print(txt)
