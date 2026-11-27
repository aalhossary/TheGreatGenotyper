#!/usr/bin/env python3
"""Compare called genotypes with truth: concordance.py truth.vcf called.vcf

Reports counts of genotypes, missing calls (./.), concordant calls, and a
confusion table of (true ALT count, called ALT count). Phase is ignored.
"""
import sys, collections
def rd(fn):
    d={}; names=None
    for l in open(fn):
        if l.startswith('##'): continue
        t=l.rstrip('\n').split('\t')
        if l[0]=='#': names=t[9:]; continue
        fmt=t[8].split(':')
        d[(t[0],int(t[1]),t[3],t[4])]=[dict(zip(fmt,x.split(':'))) for x in t[9:]]
    return names,d
tn,T=rd(sys.argv[1]); cn,C=rd(sys.argv[2])
def norm(g):
    g=g.replace('|','/')
    if '.' in g: return None
    return tuple(sorted(int(x) for x in g.split('/')))
tot=0; ok=0; miss=0; conf=collections.Counter(); gq=collections.Counter(); absent=0
for k,tv in T.items():
    if k not in C: absent+=1; continue
    for s,name in enumerate(tn):
        j=cn.index(name); a=norm(tv[s]['GT']); b=norm(C[k][j]['GT']); tot+=1
        if b is None: miss+=1; continue
        ok+= a==b; conf[(sum(a),sum(b))]+=1
print('truth sites',len(T),'called sites',len(C),'truth sites absent',absent)
print('genotypes',tot,'missing ./.',miss,'concordant',ok,'of called: %.4f'%(ok/max(1,tot-miss)),' of all: %.4f'%(ok/tot))
print('confusion (truth_altcount,called_altcount):',sorted(conf.items()))
