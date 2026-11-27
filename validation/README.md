# Validation scripts

Small-scale, end-to-end checks for The Great Genotyper on simulated data. They need no index download and run in about a minute on 4 cores.

They are a starting point for a CI test, not a finished test suite.

## What they do

| Script | Purpose |
|---|---|
| `run_small_test.sh` | Runs everything: simulate, build an index, genotype in HMM mode and in k-mer-only mode, score. Exits non-zero if concordance is below a threshold or cluster statistics are inconsistent. |
| `simulate.py` | Takes a region of `test_data/` chr21 and its variants, makes a 5-sample phased panel, and simulates samples as mosaics of panel haplotypes with 30x error-free 150 bp reads. Seeded. |
| `build_index.sh` | Builds an index by running the rules of `DatabaseBuilder/Snakefile` by hand, and records the time of each step. |
| `concordance.py` | Compares called genotypes with the simulated truth. |
| `check_cluster_stats.py` | Checks that all records of a variant cluster report the same UK and KC per sample. Regression check for mis-indexed statistics of clustered variants. |
| `kmer_counts_summary.py` | Shows whether a sequence-and-counts pair holds per-k-mer counts or one count per sequence. |
| `compare_vcfs.py` | Counts GT, GQ and GL differences between two genotyper outputs. |
| `tsan_summary.py` | Groups ThreadSanitizer reports by the first project frame. |

## Requirements

- A build of `metagraph`, `TheGreatGenotyper`, `computeHistogram` and `kmc`. Point `TGG_BUILD` at the CMake build directory (default: `../build`).
- `python3`, GNU `time` at `/usr/bin/time`, `realpath`.

Build recipe that worked on Ubuntu with system packages:

```bash
sudo apt-get install build-essential cmake gcc-12 g++-12 libjellyfish-2.0-dev jellyfish libhts-dev libbz2-dev zlib1g-dev libboost-all-dev libjemalloc-dev libdeflate-dev bcftools tabix python3-venv pkg-config time
```

```bash
CC=gcc-12 CXX=g++-12 cmake -Bbuild
```

```bash
cmake --build build -j4 --target metagraph TheGreatGenotyper PopMergeVCF computeHistogram
```

With gcc 12 the top-level `CMakeLists.txt` needs `-std=gnu++17` in place of `-std=c++17`; gcc 13 does not build the vendored libraries. Both are build problems in third-party code.

## Running

```bash
TGG_BUILD=/path/to/build bash validation/run_small_test.sh /path/to/workdir
```

Expected last lines on the default region (500 kb, 1,430 variants, 4 samples):

```
genotypes 5720 missing ./. 21 concordant 5697 of called: 0.9996  of all: 0.9960
concordance of called genotypes: 0.9996 (minimum 0.99 )
```

A larger run (5 Mb, 8 samples):

```bash
TGG_BUILD=/path/to/build bash validation/run_small_test.sh /path/to/workdir 20000000 25000000 8
```

## Deliberate differences from the committed pipeline

- KMC reads gzipped FASTA with `-fm`; the Snakefile converts to FASTQ first.
- The minimum k-mer count is 3; the committed config leaves it blank.
- The smoothing window defaults to 10,000,000 (one mean count per unitig); set `SMOOTH=1` for the Snakefile's default target.
- The genotyper is given `-e 100000000`; its default reference hash needs about 17 GB.

## Limits

The data is tiny, error-free and has a 10-haplotype panel, so the multi-subset HMM path is not exercised and timings say nothing about real scale.
