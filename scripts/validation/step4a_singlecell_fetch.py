# -*- coding: utf-8 -*-
"""Download GSE103322 (Puram 2017 HNSCC scRNA) processed matrix."""
import os, sys, ssl, urllib.request, time
sys.stdout.reconfigure(encoding='utf-8')
D=r'{DATA_ROOT}/single_cell'
os.makedirs(D, exist_ok=True)
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103322/suppl/GSE103322_HNSCC_all_data.txt.gz'
dst=os.path.join(D,'GSE103322_HNSCC_all_data.txt.gz')
ctx=ssl.create_default_context(); ctx.check_hostname=False; ctx.verify_mode=ssl.CERT_NONE
t0=time.time()
req=urllib.request.Request(url, headers={'User-Agent':'Mozilla/5.0'})
with urllib.request.urlopen(req,timeout=180,context=ctx) as r, open(dst,'wb') as f:
    n=0
    while True:
        b=r.read(1<<20)
        if not b: break
        f.write(b); n+=len(b)
print('downloaded %d bytes in %.1fs -> %s'%(n,time.time()-t0,dst))
