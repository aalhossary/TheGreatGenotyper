#!/usr/bin/env python3
"""Check per-sample statistics of clustered variants: check_cluster_stats.py out.vcf [k]

The genotyper handles variants closer than k-1 bases as one cluster and
computes one set of unique k-mers for it. Every record of a cluster must
therefore report the same UK (unique k-mers) and KC (local coverage) for a
given sample. Differences mean the statistics were read from the wrong place.
The VCF must be written with -w.

Exits 0 if all clusters are consistent, 1 otherwise, 2 if the file has no
clusters or no UK/KC fields (nothing was checked).
"""
import sys

k = int(sys.argv[2]) if len(sys.argv) > 2 else 31
clusters = []      # list of lists of (pos, per-sample dicts)
current = []
prev_chrom = None
prev_end = None
for line in open(sys.argv[1]):
    if line[0] == '#':
        continue
    t = line.rstrip('\n').split('\t')
    fmt = t[8].split(':')
    if 'UK' not in fmt or 'KC' not in fmt:
        print('no UK/KC fields; run the genotyper with -w', file=sys.stderr)
        sys.exit(2)
    start = int(t[1]) - 1
    # same rule as VariantReader: a gap of k-1 or more starts a new cluster
    if prev_chrom is not None and (t[0] != prev_chrom or start - prev_end >= k - 1):
        clusters.append(current)
        current = []
    current.append((t[0], t[1], [dict(zip(fmt, x.split(':'))) for x in t[9:]]))
    prev_chrom = t[0]
    prev_end = start + len(t[3])
if current:
    clusters.append(current)

multi = [c for c in clusters if len(c) > 1]
bad_clusters = 0
bad_records = 0
examples = []
for c in multi:
    cluster_bad = False
    for s in range(len(c[0][2])):
        first = (c[0][2][s]['UK'], c[0][2][s]['KC'])
        for chrom, pos, samples in c[1:]:
            got = (samples[s]['UK'], samples[s]['KC'])
            if got != first:
                bad_records += 1
                cluster_bad = True
                if len(examples) < 5:
                    examples.append('%s:%s sample %d: UK,KC = %s,%s but cluster starts with %s,%s at %s'
                                    % (chrom, pos, s + 1, got[0], got[1], first[0], first[1], c[0][1]))
    bad_clusters += cluster_bad

print('records', sum(map(len, clusters)), 'clusters with several variants', len(multi),
      'inconsistent clusters', bad_clusters, 'inconsistent sample-records', bad_records)
for e in examples:
    print('  ' + e)
if not multi:
    print('no clustered variants in this file; nothing was checked', file=sys.stderr)
    sys.exit(2)
sys.exit(1 if bad_clusters else 0)
