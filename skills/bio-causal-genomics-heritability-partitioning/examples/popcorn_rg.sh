#!/bin/bash
# Reference: Popcorn 1.1 (brielin/Popcorn, Python 3 port) | Verify against `popcorn --help` if
# version differs
#
# Popcorn: trans-ancestry genetic correlation from GWAS summary statistics.
# Brown 2016 AJHG 99:76.
#
# Two real bugs in the pip-installed package on this Skill's own verification environment
# (checked 2026-09-21, numpy 2.5.3 / pandas 3.0.6): `popcorn compute` crashes with
# "TypeError: only 0-dimensional arrays can be converted to Python scalars" under numpy>=2;
# fixed by installing into a venv pinned to numpy<2 (`pip install "numpy<2" "pandas<2.2" ...`).
# `popcorn fit`'s jackknife step then crashes separately with
# "AttributeError: module 'numpy' has no attribute 'bool'" (np.bool was removed in numpy>=1.24,
# and numpy<2 alone is not new enough to still have it); fixed with a one-line source patch in
# the installed package: `popcorn/jackknife.py`, `dtype=np.bool` -> `dtype=bool`.
#
# Usage:
#   bash popcorn_rg.sh <bfile1> <bfile2> <sumstats1.txt> <sumstats2.txt> <out_prefix>
#
# sumstats columns (header row): SNP A1 A2 N Z (or beta/SE, or OR/p-value -- see Popcorn's README)
# bfile1/bfile2: PLINK reference panels for population 1 and 2 (real cross-ancestry panels;
# same-ancestry panels make the model's h2 estimates unstable/inflated -- see Common Errors below)

set -euo pipefail

BFILE1=${1:?Usage: popcorn_rg.sh <bfile1> <bfile2> <sumstats1.txt> <sumstats2.txt> <out_prefix>}
BFILE2=${2:?provide bfile2}
SUMSTATS1=${3:?provide sumstats1.txt}
SUMSTATS2=${4:?provide sumstats2.txt}
OUT=${5:?provide out_prefix}

# Step 1: compute cross-population LD/covariance scores from the two reference panels
popcorn compute -v 1 --bfile1 "${BFILE1}" --bfile2 "${BFILE2}" "${OUT}_scores.txt"

# Step 2: fit heritability + trans-ancestry genetic correlation to the summary statistics
popcorn fit -v 1 --cfile "${OUT}_scores.txt" \
    --sfile1 "${SUMSTATS1}" --sfile2 "${SUMSTATS2}" \
    "${OUT}_result.txt"

echo "Result in ${OUT}_result.txt: h1^2, h2^2 (per-population observed-scale heritability),"
echo "pgi (genetic impact correlation), each with SE, Z, and P."
echo ""
echo "Common error: h1^2/h2^2 wildly outside [0,1] (verified 2026-09-21 on a same-ancestry toy"
echo "panel split in two -- gave h2 ~ 150) means the two bfiles do not have genuinely distinct"
echo "population-level LD/allele-frequency structure; Popcorn's cross-covariance-score model"
echo "assumes real cross-ancestry panels (e.g. 1000G EUR vs EAS), not two subsets of one cohort."
