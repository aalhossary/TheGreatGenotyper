#!/bin/bash
# Build a small index by running the rules of DatabaseBuilder/Snakefile by hand.
#
# usage: build_index.sh <simdir> <workdir> <smoothing_window> "<sample names>"
#   <simdir>   directory with <sample>.fa.gz read files (see simulate.py)
#   <workdir>  output directory; the index ends up in <workdir>/idx/
#   <smoothing_window>  1 = per-k-mer counts (the Snakefile default target);
#                       10000000 = one mean count per unitig
# environment:
#   TGG_BUILD  CMake build directory (default: <repo>/build)
#   THREADS    threads per step (default: 4)
#   MEM_GB     --mem-cap-gb for metagraph (default: 4)
#
# Differences from the Snakefile, on purpose:
#   - KMC reads gzipped FASTA directly with -fm (no seqtk conversion to FASTQ)
#   - kmcMinCount is fixed at 3 (the committed config leaves it blank)
#   - the Snakefile passes the smoothing value as --subsample; so does this
# Step timings are written to <workdir>/timings.tsv.
set -uo pipefail
SIM=$(realpath "$1"); D=$2; sm=$3; SAMPLES=$4
HERE=$(cd "$(dirname "$0")" && pwd)
B=$(realpath "${TGG_BUILD:-$HERE/../build}"); P=${THREADS:-4}; MEM=${MEM_GB:-4}
MG=$B/metagraph; KMC=$B/KMC/kmc; HISTO=$B/DatabaseBuilder/computeHistogram
for x in "$MG" "$KMC" "$HISTO"; do [ -x "$x" ] || { echo "missing binary: $x" >&2; exit 1; }; done
mkdir -p "$D"; cd "$D"
T(){ local name=$1; shift; /usr/bin/time -o .time.tmp -f "%e\t%U\t%S\t%M" "$@"; local rc=$?; printf "%s\trc=%s\t%s\n" "$name" "$rc" "$(tail -1 .time.tmp)" >> timings.tsv; return $rc; }
echo -e "step\trc\twall_s\tuser_s\tsys_s\tmaxrssKB" > timings.tsv
mkdir -p result/kmc result/histograms result/graphs result/unitigs/smooth_$sm tmp
for s in $SAMPLES; do
  mkdir -p tmp/kmc_$s; echo $SIM/$s.fa.gz > result/kmc/$s.lst
  T kmc_$s $KMC -ci3 -t$P -k31 -m16 -fm @result/kmc/$s.lst result/kmc/$s tmp/kmc_$s/ &> result/kmc/$s.log
  T histo_$s $HISTO -i result/kmc/$s -m 1000 -o result/histograms/$s.histo &>> result/kmc/$s.log
  rm -rf tmp/kmc_$s; mkdir -p tmp/graph_$s; cp result/kmc/$s.kmc_pre result/kmc/$s.kmc_suf tmp/graph_$s/
  T kmc_to_graph_$s $MG build -k 31 --mode canonical --count-kmers --count-width 16 --mem-cap-gb $MEM --disk-swap tmp/graph_$s/ -p $P -o tmp/graph_$s/graph tmp/graph_$s/$s.kmc_suf &> result/graphs/$s.log
  mv tmp/graph_$s/graph.dbg result/graphs/$s.dbg; mv tmp/graph_$s/graph.dbg.weights result/graphs/$s.dbg.weights; rm -rf tmp/graph_$s
  mkdir -p tmp/clean_$s; cp result/graphs/$s.dbg result/graphs/$s.dbg.weights tmp/clean_$s/
  T clean_$s $MG clean --to-fasta --primary-kmers --prune-tips 0 --prune-unitigs 1 --smoothing-window $sm -p $P -o tmp/clean_$s/$s tmp/clean_$s/$s.dbg &> result/unitigs/smooth_$sm/$s.log
  mv tmp/clean_$s/$s.fasta.gz tmp/clean_$s/$s.kmer_counts.gz result/unitigs/smooth_$sm/; rm -rf tmp/clean_$s
done
out=result/smooth_$sm; mkdir -p $out tmp/main
U=$(for s in $SAMPLES; do echo result/unitigs/smooth_$sm/$s.fasta.gz; done)
L=$out/buildMainGraphlog
T main_build_canonical $MG build -k 31 --mode canonical --mem-cap-gb $MEM --disk-swap tmp/main/ -p $P -o tmp/main/graph_canonical $U &> $L
T main_transform_primary $MG transform --to-fasta --primary-kmers -p $P -o tmp/main/graph_primary tmp/main/graph_canonical.dbg &>> $L
rm tmp/main/graph_canonical.dbg
T main_build_primary $MG build -k 31 --mode primary --mem-cap-gb $MEM --disk-swap tmp/main/ -p $P -o $out/graph.dbg tmp/main/graph_primary.fasta.gz &>> $L
rm -rf tmp/main
mkdir -p $out/columns
T annotate $MG annotate -i $out/graph.dbg --anno-filename --separately --count-kmers --count-width 16 -o $out/columns/ --threads-each 8 -p 4 $U &> $out/annotate.log
C=$(for s in $SAMPLES; do echo $out/columns/$s.fasta.gz.column.annodbg; done)
mkdir -p tmp/opt $out/rowDiff; L=$out/optimizeAnnotationColumns.log; : > $L
for st in 0 1 2; do
  T rowdiff_stage$st $MG transform_anno --anno-type row_diff --count-kmers --row-diff-stage $st --mem-cap-gb $MEM --disk-swap tmp/opt -i $out/graph.dbg -o $out/rowDiff/out -p $P $C &>> $L
done
find $out/rowDiff/ -name "*.column.annodbg" > $out/rowdiff_columns.txt
T brwt bash -c "cat $out/rowdiff_columns.txt | $MG transform_anno --anno-type row_diff_int_brwt --greedy --fast --subsample $sm -i $out/graph.dbg -o $out/annotation -p $P --parallel-nodes 10" &>> $L
T relax $MG relax_brwt -p $P --relax-arity 32 -o $out/annotation.relaxed $out/annotation.row_diff_int_brwt.annodbg &>> $L
# graph.desc.tsv: sample name, k-mer coverage (from the histogram), annotation label
: > $out/graph.desc.tsv
for s in $SAMPLES; do
  printf "%s\t%s\t%s\n" "$s" "$(grep parameters result/histograms/$s.histo | cut -f3)" "result/unitigs/smooth_$sm/$s.fasta.gz" >> $out/graph.desc.tsv
done
mkdir -p idx; cp $out/graph.dbg $out/graph.desc.tsv $out/annotation.relaxed.row_diff_int_brwt.annodbg idx/
echo "$PWD/idx/" > indexes
[ -s idx/graph.dbg ] && [ -s idx/annotation.relaxed.row_diff_int_brwt.annodbg ] && [ -s idx/graph.desc.tsv ] || { echo "index build failed; see logs under $PWD/result" >&2; exit 1; }
echo "index written to $PWD/idx/"
