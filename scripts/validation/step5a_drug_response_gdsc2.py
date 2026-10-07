# -*- coding: utf-8 -*-
"""External validation #3: drug connectivity.
Does the inflammatory-core expression score predict sensitivity to drugs across
GDSC2 cell lines (and specifically in head-and-neck lines)?"""
import os, sys, re, csv
import numpy as np, pandas as pd
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')
DEP=r'{DATA_ROOT}/cancer_application/data/depmap'
out=[]

# --- DepMap expression: core score ---
EXPR=os.path.join(DEP,'OmicsExpressionProteinCodingGenesTPMLogp1.csv')
hdr=open(EXPR,encoding='utf-8',errors='ignore').readline().rstrip('\n').split(',')
CORE=['STAT1','JUN','FOSL1','RELA','CXCL8','CXCL10','MMP9','PTGS2','E2F1','ETS1']
def find(g):
    for i,c in enumerate(hdr):
        if c.startswith(g+' ('): return i
    return None
idx={g:find(g) for g in CORE}
idx={g:i for g,i in idx.items() if i is not None}
out.append('core genes in DepMap expression: %d'%len(idx))
cols=[0]+list(idx.values())
ex=pd.read_csv(EXPR, usecols=cols, low_memory=False)
ex.columns=['ModelID']+list(idx.keys())
ex=ex.set_index('ModelID')
z=((ex-ex.mean())/ex.std())
ex['CORE_SCORE']=z.mean(axis=1)
out.append('expression matrix shape=%s'%str(ex.shape))

model=pd.read_csv(os.path.join(DEP,'Model.csv'), low_memory=False)
model['key']=model['CellLineName'].astype(str).str.upper().str.replace(r'[^A-Z0-9]','',regex=True)
mod=model.set_index('ModelID')
ex['lineage']=[mod['OncotreeLineage'].get(i,'') for i in ex.index]
ex['cellname']=[mod['CellLineName'].get(i,'') for i in ex.index]
ex['key']=[mod['key'].get(i,'') for i in ex.index]

# --- GDSC2 ---
gd=pd.read_csv(os.path.join(DEP.replace('depmap','gdsc'),'GDSC2-dataset.csv'), low_memory=False)
out.append('GDSC2 shape=%s'%str(gd.shape))
gd['key']=gd['CELL_LINE_NAME'].astype(str).str.upper().str.replace(r'[^A-Z0-9]','',regex=True)
mg=gd.merge(ex[['key','CORE_SCORE','lineage']], on='key', how='inner')
out.append('GDSC2 x DepMap matched rows=%d  unique lines=%d  drugs=%d'%(
    len(mg), mg['key'].nunique(), mg['DRUG_NAME'].nunique()))

def corr_table(df, minn=20):
    rows=[]
    for drug,d in df.groupby('DRUG_NAME'):
        d=d.dropna(subset=['CORE_SCORE','AUC'])
        if len(d)<minn: continue
        rho,p=stats.spearmanr(d['CORE_SCORE'], d['AUC'])
        rows.append((drug, d['PATHWAY_NAME'].iloc[0], d['PUTATIVE_TARGET'].iloc[0], len(d), rho, p))
    t=pd.DataFrame(rows, columns=['drug','pathway','target','n','rho','p']).sort_values('rho')
    return t

t_all=corr_table(mg)
t_all['q']=stats.false_discovery_control(t_all['p'].values)
out.append('\n=== ALL lines: core-score vs AUC (rho<0 => core-high = more sensitive) ===')
out.append('drugs tested=%d'%len(t_all))
out.append('-- most SENSITIVE (most negative rho) --')
for _,r in t_all.head(15).iterrows():
    out.append('  %-28s rho=%+.3f p=%.1e q=%.2e  n=%d  [%s]'%(r['drug'],r['rho'],r['p'],r['q'],r['n'],str(r['target'])[:18]))
out.append('-- most RESISTANT (most positive rho) --')
for _,r in t_all.tail(10).iterrows():
    out.append('  %-28s rho=%+.3f p=%.1e q=%.2e  n=%d  [%s]'%(r['drug'],r['rho'],r['p'],r['q'],r['n'],str(r['target'])[:18]))

# highlight inflammation-pathway drugs
mask=t_all['pathway'].astype(str).str.contains('NF|JAK|STAT|Apoptosis|Immune|Inflammation|TNF', case=False, na=False)
out.append('\n-- inflammation/apoptosis-pathway drugs --')
for _,r in t_all[mask].iterrows():
    out.append('  %-28s rho=%+.3f p=%.1e  n=%d  pw=%s'%(r['drug'],r['rho'],r['p'],r['n'],str(r['pathway'])[:22]))

# head&neck subset
hn=mg[mg['lineage'].astype(str).str.contains('Head and Neck',case=False,na=False)]
out.append('\n=== HEAD&NECK lines only ===  rows=%d lines=%d'%(len(hn), hn['key'].nunique()))
t_hn=corr_table(hn, minn=8)
out.append('drugs tested=%d'%len(t_hn))
for _,r in t_hn.head(12).iterrows():
    out.append('  %-28s rho=%+.3f p=%.1e  n=%d  [%s]'%(r['drug'],r['rho'],r['p'],r['n'],str(r['target'])[:18]))

# specific candidate drugs
out.append('\n=== named candidate drugs (all lines) ===')
for name in ['Bortezomib','Cisplatin','5-Fluorouracil','Docetaxel','Methotrexate','Nutlin','Ruxolitinib','Ibrutinib','Dexamethasone','Vorinostat']:
    hit=t_all[t_all['drug'].astype(str).str.contains(name, case=False, na=False)]
    for _,r in hit.iterrows():
        out.append('  %-28s rho=%+.3f p=%.1e n=%d'%(r['drug'],r['rho'],r['p'],r['n']))

open(r'./_gdsc.txt','w',encoding='utf-8').write('\n'.join(out))
print('\n'.join(out)[:5000])
