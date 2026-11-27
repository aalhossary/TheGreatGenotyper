#!/bin/bash
# End-to-end check on simulated data: simulate, build an index, genotype, score.
#
# usage: run_small_test.sh <workdir> [start] [end] [n_samples]
#   default region: chr21:20,000,000-20,500,000 of test_data, 4 samples
# environment: TGG_BUILD (default <repo>/build), THREADS (default 4),
#              SMOOTH (default 10000000), MIN_CONCORDANCE (default 0.99)
# Exits non-zero if concordance of called genotypes is below MIN_CONCORDANCE.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/.." && pwd)
D=$1; START=${2:-20000000}; END=${3:-20500000}; NS=${4:-4}
B=$(realpath "${TGG_BUILD:-$REPO/build}"); P=${THREADS:-4}; SM=${SMOOTH:-10000000}
GG=$B/pangenie/src/TheGreatGenotyper
[ -x "$GG" ] || { echo "missing binary: $GG" >&2; exit 1; }
mkdir -p "$D/sim"; D=$(realpath "$D")
python3 "$HERE/simulate.py" "$REPO/test_data/GRCh38_chr21.fa.gz" "$REPO/test_data/test.vcf.gz" "$D/sim" "$START" "$END" "$NS"
SAMPLES=$(for i in $(seq 1 "$NS"); do printf "S%d " "$i"; done)
TGG_BUILD=$B THREADS=$P bash "$HERE/build_index.sh" "$D/sim" "$D/index" "$SM" "$SAMPLES"
mkdir -p "$D/run"; cd "$D/run"
# -e sets the reference hash size; the default needs about 17 GB.
"$GG" -g -i "$D/index/indexes" -j "$P" -t "$P" -e 100000000 -r "$D/sim/ref.fa" -y "$D/run/emissions" -v "$D/sim/panel.vcf" -o "$D/run/out.vcf" 2> "$D/run/log"
python3 "$HERE/concordance.py" "$D/sim/truth.vcf" "$D/run/out.vcf" | tee "$D/run/concordance.txt"
c=$(sed -n 's/.*of called: \([0-9.]*\).*/\1/p' "$D/run/concordance.txt")
python3 -c "import sys; c=float('$c'); m=float('${MIN_CONCORDANCE:-0.99}'); print('concordance of called genotypes:', c, '(minimum', m, ')'); sys.exit(0 if c>=m else 1)"
