#!/usr/bin/env python3
"""Summarise a MetaGraph sequence + counts pair: kmer_counts_summary.py <prefix>

Reads <prefix>.fasta.gz and <prefix>.kmer_counts.gz (4-byte k, then one uint32
per k-mer) and reports how many sequences carry a single constant count.
Used to see what --smoothing-window does.
"""
import gzip, struct, sys, collections
pre=sys.argv[1]
c=gzip.open(pre+'.kmer_counts.gz','rb').read()
k=struct.unpack('<I',c[:4])[0]; off=4; nseq=0; nconst=0; nk=0; vals=collections.Counter(); ex=None; multi=0
for l in gzip.open(pre+'.fasta.gz','rt'):
    if l[0]=='>': continue
    n=len(l.strip())-k+1; v=struct.unpack('<%dI'%n,c[off:off+4*n]); off+=4*n
    nseq+=1; nk+=n
    if n>1: multi+=1
    if len(set(v))==1: nconst+=1
    elif ex is None and n>40: ex=v[:25]
    vals.update(v)
print(pre,'k',k,'seqs',nseq,'kmers',nk,'seqs_len>1kmer',multi,'seqs_constant_count',nconst,'distinct_count_values',len(vals),'bytes_left',len(c)-off)
print('  example non-constant:',ex); print('  most common values',vals.most_common(6))
