# -*- coding: utf-8 -*-
"""GSE103322 (Puram 2017 HNSCC scRNA): is the inflammatory core expressed in
malignant cells, or confined to immune/stromal cells?
ROW3 'classified  as cancer cell' (1=malignant), ROW5 'non-cancer cell type'."""
import os, sys, gzip, collections
import numpy as np
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')
P=r'{DATA_ROOT}/single_cell/GSE103322_HNSCC_all_data.txt.gz'
out=[]

CORE=['STAT1','JUN','FOSL1','RELA','CXCL8','CXCL10','MMP9','PTGS2','E2F1','ETS1']
MARK=['EPCAM','PTPRC','CD3D','CD68','COL1A1','MS4A1','KRT14','KRT5','LYZ','PECAM1']
WANT=set(CORE+MARK)
ANNOT_KEEP={'classified  as cancer cell','non-cancer cell type','classified as non-cancer cells'}

header=None; gene={}; annot={}
with gzip.open(P,'rt',encoding='utf-8',errors='ignore') as f:
    for line in f:
        if header is None:
            header=[c.strip().strip('"') for c in line.rstrip('\n').split('\t')]; continue
        i=line.find('\t'); lab=line[:i].strip().strip('"').strip("'").strip()
        if lab in WANT:
            v=[x.strip().strip("'").strip('"') for x in line[i+1:].rstrip('\n').split('\t')]
            gene[lab]=np.array([float(x) if x not in('','NA') else np.nan for x in v],dtype=np.float32)
        elif lab in ANNOT_KEEP:
            annot[lab]=[x.strip().strip("'").strip('"') for x in line[i+1:].rstrip('\n').split('\t')]

maln=np.array(annot['classified  as cancer cell'])
ctyr=annot['non-cancer cell type']
ctype=np.array(['Malignant' if m=='1' else (ctyr[i] if ctyr[i] not in ('0','') else 'Unassigned')
                for i,m in enumerate(maln)])
nc=collections.Counter(ctype)
out.append('cells=%d  cell types=%s'%(len(ctype),dict(nc)))

core=[g for g in CORE if g in gene]
Z=np.vstack([(gene[g]-np.nanmean(gene[g]))/(np.nanstd(gene[g])+1e-9) for g in core])
score=np.nanmean(Z,axis=0)
out.append('core genes used=%s (of %d)'%(core,len(CORE)))

out.append('\n=== mean core z-score by cell type ===')
rows=[]
for c,n in nc.items():
    idx=np.where(ctype==c)[0]
    rows.append((c,n,float(np.nanmean(score[idx])),float(np.nanmedian(score[idx]))))
rows.sort(key=lambda r:-r[2])
out.append('%-14s %6s %10s %10s'%('celltype','n','mean_z','median_z'))
for c,n,m,md in rows:
    out.append('%-14s %6d %10.3f %10.3f'%(c,n,m,md))

mi=np.where(ctype=='Malignant')[0]; ti=np.where(ctype=='T cell')[0]
fi=np.where(ctype=='Fibroblast')[0]; ma=np.where(ctype=='Macrophage')[0]
out.append('\nMalignant n=%d  mean core z=%+.3f'%(len(mi),np.nanmean(score[mi])))
out.append('MWU Malignant vs T cell      p=%.2e'%stats.mannwhitneyu(score[mi],score[ti]).pvalue)
out.append('MWU Malignant vs Fibroblast  p=%.2e'%stats.mannwhitneyu(score[mi],score[fi]).pvalue)
out.append('MWU Malignant vs Macrophage  p=%.2e'%stats.mannwhitneyu(score[mi],score[ma]).pvalue)

# per-gene detection and mean level, malignant vs T cell vs fibroblast
out.append('\n=== per-gene: fraction expressing (>0) and mean level, by cell type ===')
out.append('%-8s | %-19s | %-19s | %-19s'%('gene','Malignant','T cell','Fibroblast'))
out.append('%-8s | %8s %10s | %8s %10s | %8s %10s'%('','frac>0','mean','frac>0','mean','frac>0','mean'))
for g in core:
    v=gene[g]
    def st(idx): return float(np.nanmean(v[idx]>0)), float(np.nanmean(v[idx]))
    fm,mm=st(mi); ft,mt=st(ti); ff,mf=st(fi)
    out.append('%-8s | %8.3f %10.2f | %8.3f %10.2f | %8.3f %10.2f'%(g,fm,mm,ft,mt,ff,mf))

# malignant cells: fraction with core score above global median
med=np.nanmedian(score)
out.append('\nfraction of cells with core_z above global median:')
for c,n in nc.most_common(8):
    idx=np.where(ctype==c)[0]
    out.append('  %-14s %.3f'%(c,float(np.nanmean(score[idx]>med))))

open(r'./_sc.txt','w',encoding='utf-8').write('\n'.join(out))
print('\n'.join(out))
