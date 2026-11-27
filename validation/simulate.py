#!/usr/bin/env python3
"""Simulate a small phased panel and short-read samples from a reference region.

Writes into <outdir>: ref.fa, panel.vcf (phased panel, samples P1..), truth.vcf
(true genotypes of simulated samples S1..), and S<n>.fa.gz (error-free 150 bp
reads, about 30x, single-line FASTA). Sample haplotypes are mosaics of panel
haplotypes. Only biallelic, non-overlapping variants are kept. Seeded, so
output is reproducible.
"""
# usage: simulate.py chr21.fa.gz test.vcf.gz outdir start end
import sys, gzip, random
fa, vcf, out, start, end = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
random.seed(42)
seq=[]
with gzip.open(fa,'rt') as f:
    for l in f:
        if l[0]=='>': continue
        seq.append(l.strip())
ref=''.join(seq)[start:end].upper()
assert 'N' not in ref, ref.count('N')
with open(out+'/ref.fa','w') as f:
    f.write('>chr21\n')
    for i in range(0,len(ref),60): f.write(ref[i:i+60]+'\n')
vars_=[]; last_end=0
with gzip.open(vcf,'rt') as f:
    for l in f:
        if l[0]=='#': continue
        t=l.split('\t'); pos=int(t[1])-1
        if pos<start+200 or pos>end-200: continue
        if ',' in t[4]: continue
        p=pos-start; r=t[3].upper(); a=t[4].upper()
        if p+len(r)>len(ref)-200: continue
        if ref[p:p+len(r)]!=r: continue
        if p<=last_end+1: continue   # non-overlapping, non-adjacent
        if set(a)-set('ACGT'): continue
        vars_.append((p,r,a)); last_end=p+len(r)
NP=5   # panel samples -> 10 haplotypes
panel=[]   # per variant list of 10 alleles
for v in vars_:
    af=random.choice([0.1,0.2,0.3,0.5,0.7])
    h=[1 if random.random()<af else 0 for _ in range(2*NP)]
    if sum(h)==0: h[random.randrange(2*NP)]=1
    panel.append(h)
NS=int(sys.argv[6]) if len(sys.argv)>6 else 4
truth=[]  # per sample: 2 haplotype allele lists
for s in range(NS):
    haps=[]
    for hh in range(2):
        src=random.randrange(2*NP); al=[]; nextsw=random.expovariate(1/50000.)
        for i,v in enumerate(vars_):
            while v[0]>nextsw:
                src=random.randrange(2*NP); nextsw+=random.expovariate(1/50000.)
            al.append(panel[i][src])
        haps.append(al)
    truth.append(haps)
hdr='##fileformat=VCFv4.2\n##contig=<ID=chr21,length=%d>\n##FILTER=<ID=PASS,Description="All filters passed">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'%len(ref)
with open(out+'/panel.vcf','w') as f:
    f.write(hdr+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join('P%d'%i for i in range(1,NP+1))+'\n')
    for v,h in zip(vars_,panel):
        f.write('chr21\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT\t'%(v[0]+1,v[1],v[2])+'\t'.join('%d|%d'%(h[2*i],h[2*i+1]) for i in range(NP))+'\n')
with open(out+'/truth.vcf','w') as f:
    f.write(hdr+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join('S%d'%i for i in range(1,NS+1))+'\n')
    for i,v in enumerate(vars_):
        f.write('chr21\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT\t'%(v[0]+1,v[1],v[2])+'\t'.join('%d|%d'%(truth[s][0][i],truth[s][1][i]) for s in range(NS))+'\n')
comp=str.maketrans('ACGT','TGCA')
RL=150; COV=15
for s in range(NS):
    with gzip.open(out+'/S%d.fa.gz'%(s+1),'wt',compresslevel=1) as f:
        n=0
        for hh in range(2):
            parts=[]; prev=0
            for i,v in enumerate(vars_):
                parts.append(ref[prev:v[0]])
                parts.append(v[2] if truth[s][hh][i] else v[1]); prev=v[0]+len(v[1])
            parts.append(ref[prev:]); hs=''.join(parts)
            for _ in range(len(hs)*COV//RL):
                p=random.randrange(len(hs)-RL+1); r=hs[p:p+RL]
                if random.random()<0.5: r=r.translate(comp)[::-1]
                n+=1; f.write('>r%d\n%s\n'%(n,r))
print('variants',len(vars_),'reflen',len(ref))
