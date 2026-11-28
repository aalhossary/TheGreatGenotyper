#!/bin/bash
# Regression test for PopMergeVCF: the GQ threshold is computed per variant, so
# the output for a line must be the same whether the line is processed alone
# or after other lines. It used to depend on all earlier lines, because the
# per-line quality lists were cleared on a copy.
#
# usage: test_line_independence.sh <path to PopMergeVCF binary>
set -euo pipefail
BIN=$1
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT

# Six samples. Line 1 has high qualities, line 2 low ones, line 3 a genotype
# class (1/1) that the earlier lines do not have, line 4 lacks classes that
# earlier lines had.
{
printf '##fileformat=VCFv4.2\n'
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6\n'
printf 'chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT:GQ\t0/0:900\t0/0:800\t0/0:700\t0/0:600\t0/1:500\t0/1:400\n'
printf 'chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT:GQ\t0/0:90\t0/0:80\t0/0:70\t0/0:60\t0/1:50\t0/1:40\n'
printf 'chr1\t300\t.\tT\tA\t.\tPASS\t.\tGT:GQ\t1/1:30\t1/1:20\t1/1:10\t1/1:5\t0/1:45\t0/1:35\n'
printf 'chr1\t400\t.\tC\tG\t.\tPASS\t.\tGT:GQ\t1/1:300\t1/1:200\t1/1:100\t1/1:50\t1/1:25\t1/1:12\n'
} > "$T/in.vcf"

"$BIN" < "$T/in.vcf" > "$T/all.out"

# the same lines, each processed on its own
grep '^#' "$T/in.vcf" > "$T/single.out"
grep -v '^#' "$T/in.vcf" | while IFS= read -r line; do
    printf '%s\n' "$line" | "$BIN" >> "$T/single.out"
done

if ! diff "$T/all.out" "$T/single.out" > "$T/diff"; then
    echo "FAIL: output of a line depends on earlier lines" >&2
    cat "$T/diff" >&2
    exit 1
fi

# expected thresholds: the value at position size/4 of each genotype class
expected_line2=$'chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT:GQ\t0/0:90\t0/0:80\t0/0:70\t./.:60\t0/1:50\t0/1:40'
if [ "$(sed -n 4p "$T/all.out")" != "$expected_line2" ]; then
    echo "FAIL: unexpected output for the second variant" >&2
    sed -n 4p "$T/all.out" >&2
    exit 1
fi
echo "PopMergeVCF line independence: passed"
