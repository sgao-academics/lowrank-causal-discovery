# -*- coding: utf-8 -*-
"""External validation #4: is the 'inflammatory core' just immune/stromal infiltration?
Marker-based deconvolution of the two Indian cohorts + adjusted test."""
import os, sys, gzip
import numpy as np
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')
OUT=r'{DATA_ROOT}/cancer_application/data/indian_oral'
ANNOT=os.path.join(OUT,'GPL6480.annot.gz')
out=[]

# --- platform annotation ---
probe2sym={}; hdr=None
with gzip.open(ANNOT,'rt',encoding='utf-8',errors='ignore') as f:
    for line in f:
        if not line or line[0] in '#!^': continue
        cols=[c.strip().strip('"') for c in line.rstrip('\n').split('\t')]
        if hdr is None:
            if cols and cols[0]=='ID':
                hdr=cols
                for cand in ('GENE_SYMBOL','Gene symbol','GENE','symbol','GENE_NAME'):
                    if cand in hdr: isym=hdr.index(cand); break
                ip=hdr.index('ID')
            continue
        if len(cols)>max(ip,isym):
            s=cols[isym].strip()
            if s: probe2sym[cols[ip].strip()]=s.upper()
out.append('probes with symbol: %d'%len(probe2sym))

GROUPS={'GSE85195':('Oral Tumor','Leukoplakia'),'GSE23558':('Oral Tumor','Normal')}
IMM=['PTPRC','CD3D','CD3E','CD2','CD8A','CD4','CD19','MS4A1','CD79A','CD79B','IGHM','CD68','CD14','CD163','CSF1R','ITGAM','FCGR3A','NKG7','GNLY','PRF1','IL2RA','CD27','CXCR4']
STR=['COL1A1','COL1A2','COL3A1','COL5A1','COL5A2','COL6A1','COL6A2','COL6A3','FAP','PDGFRB','DCN','LUM','ACTA2','TAGLN','POSTN','THY1','FN1','MMP2','VIM','SPARC','COL1A1']
CORE=['STAT1','JUN','FOSL1','RELA','CXCL8','CXCL10','MMP9','PTGS2','E2F1','ETS1']

def load(acc,case_kw,ctrl_kw):
    p=os.path.join(OUT,acc+'_series_matrix.txt.gz')
    sources=[];val=[];pids=[];in_tbl=False
    with gzip.open(p,'rt',encoding='utf-8',errors='ignore') as f:
        for line in f:
            if line.startswith('!Sample_source_name'):
                sources=[x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]]
            elif line.startswith('!series_matrix_table_begin'):
                in_tbl=True; f.readline(); continue
            elif line.startswith('!series_matrix_table_end'):
                in_tbl=False
            elif in_tbl:
                c=line.rstrip('\n').split('\t')
                pids.append(c[0].strip('"')); val.append(c[1:])
    M=np.array([[float(v) if v not in('','NA') else np.nan for v in r] for r in val],dtype=np.float32)
    gidx={}
    for i,pr in enumerate(pids):
        s=probe2sym.get(pr)
        if s: gidx.setdefault(s,[]).append(i)
    genes=[];vecs=[]
    for g,ii in gidx.items():
        sub=M[ii]; genes.append(g); vecs.append(sub[int(np.nanargmax(np.nanmean(sub,axis=1)))])
    E=np.array(vecs,dtype=np.float32)
    gi={g:i for i,g in enumerate(genes)}
    lab=np.array(['case' if case_kw in s else ('ctrl' if ctrl_kw in s else 'other') for s in sources])
    return E,gi,lab,sources

def score(E,gi,geneset):
    pres=[g for g in geneset if g in gi]
    z=[]
    for g in pres:
        v=E[gi[g]].astype(float); v=(v-np.nanmean(v))/(np.nanstd(v)+1e-9); z.append(v)
    return np.nanmean(np.vstack(z),axis=0), pres

def zscore_rows(E,gi,genes):
    return np.vstack([ (E[gi[g]]-np.nanmean(E[gi[g]]))/(np.nanstd(E[gi[g]])+1e-9) for g in genes if g in gi ])

for acc,(ck,nk) in GROUPS.items():
    E,gi,lab,src=load(acc,ck,nk)
    ci=np.where(lab=='case')[0]; ni=np.where(lab=='ctrl')[0]
    out.append('\n########## %s : %d case vs %d ctrl (genes=%d) ##########'%(acc,len(ci),len(ni),len(gi)))
    core,gc=score(E,gi,CORE); imm,im=score(E,gi,IMM); st,sp=score(E,gi,STR)
    out.append('core genes used=%s'%gc)
    out.append('immune markers=%d  stromal markers=%d'%(len(im),len(sp)))
    # confound exists?
    out.append('immune score: case %.3f vs ctrl %.3f  p=%.2e'%(np.nanmean(imm[ci]),np.nanmean(imm[ni]),stats.mannwhitneyu(imm[ci],imm[ni]).pvalue))
    out.append('stromal score: case %.3f vs ctrl %.3f  p=%.2e'%(np.nanmean(st[ci]),np.nanmean(st[ni]),stats.mannwhitneyu(st[ci],st[ni]).pvalue))
    out.append('core  score: case %.3f vs ctrl %.3f  p=%.2e'%(np.nanmean(core[ci]),np.nanmean(core[ni]),stats.mannwhitneyu(core[ci],core[ni]).pvalue))
    r_i=stats.spearmanr(core,imm).statistic; r_s=stats.spearmanr(core,st).statistic
    out.append('Spearman core~immune = %.3f ; core~stroma = %.3f'%(r_i,r_s))
    # adjusted: residualise core on immune+stroma (linear), then case vs ctrl
    X=np.column_stack([np.ones_like(imm),imm,st])
    ok=np.isfinite(core)&np.isfinite(imm)&np.isfinite(st)
    beta,_,_,_=np.linalg.lstsq(X[ok],core[ok],rcond=None)
    resid=np.full_like(core,np.nan); resid[ok]=core[ok]-X[ok]@beta
    out.append('ADJUSTED core (regressed on immune+stroma): case %.3f vs ctrl %.3f  p=%.2e'%(
        np.nanmean(resid[ci]),np.nanmean(resid[ni]),stats.mannwhitneyu(resid[ci],resid[ni]).pvalue))
    # per-gene adjusted
    out.append('per-gene  raw p / adjusted p (group coef after immune+stroma):')
    for g in CORE:
        if g not in gi: continue
        v=E[gi[g]].astype(float)
        ok=np.isfinite(v)&np.isfinite(imm)&np.isfinite(st)
        grp=np.zeros(len(v)); grp[ci]=1
        p_raw=stats.mannwhitneyu(v[ci],v[ni]).pvalue
        Xd=np.column_stack([np.ones(ok.sum()),grp[ok],imm[ok],st[ok]])
        b,_,_,_=np.linalg.lstsq(Xd,v[ok],rcond=None)
        yhat=Xd@b; r=v[ok]-yhat
        dof=ok.sum()-Xd.shape[1]
        s2=(r@r)/dof
        XtXinv=np.linalg.inv(Xd.T@Xd)
        se=np.sqrt(s2*XtXinv[1,1]); t=b[1]/se
        p_adj=2*stats.t.sf(abs(t),dof)
        out.append('   %-7s raw=%.2e  adj=%.2e  (coef=%+.2f)'%(g,p_raw,p_adj,b[1]))

open(r'./_deconv.txt','w',encoding='utf-8').write('\n'.join(out))
print('\n'.join(out))
