#!/usr/bin/env python3
"""Compare two genotyper outputs field by field: compare_vcfs.py a.vcf b.vcf label

Counts differing GT and GQ values and the largest GL difference. Both files
must hold the same records in the same order.
"""
import sys
def rd(fn):
    out=[]
    for l in open(fn):
        if l[0]=='#': continue
        t=l.rstrip('\n').split('\t'); f=t[8].split(':')
        out.append([dict(zip(f,x.split(':'))) for x in t[9:]])
    return out
a=rd(sys.argv[1]); b=rd(sys.argv[2]); assert len(a)==len(b)
n=gt=0; mx=0; gqd=0; glmax=0.0; big=0; inf=0
for ra,rb in zip(a,b):
    for x,y in zip(ra,rb):
        n+=1; gt+= x['GT']!=y['GT']
        d=abs(int(x['GQ'])-int(y['GQ'])); gqd+= d>0; mx=max(mx,d)
        if 'GL' in x:
            for p,q in zip(x['GL'].split(','),y['GL'].split(',')):
                try:
                    p=float(p); q=float(q)
                    if p!=p or q!=q or abs(p)==float('inf') or abs(q)==float('inf'):
                        inf+= (str(p)!=str(q)); continue
                    glmax=max(glmax,abs(p-q))
                except: pass
print('%s: genotypes=%d differing_GT=%d differing_GQ=%d max_GQ_diff=%d max_abs_GL_diff=%.4g inf/nan_mismatch=%d'%(sys.argv[3],n,gt,gqd,mx,glmax,inf))
