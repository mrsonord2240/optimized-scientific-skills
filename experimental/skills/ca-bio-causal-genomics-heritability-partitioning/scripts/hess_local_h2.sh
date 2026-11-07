#!/bin/bash
# Reference: HESS 0.5.4-beta (huwenboshi/hess) | Verify against `hess.py --help` if version differs
#
# HESS local (per-locus) SNP-heritability, step 1 (run once per chromosome).
# Shi 2016 AJHG 99:139: per-locus h2 via quadratic form on LD-projected effect estimates.
#
# The shipped huwenboshi/hess code does NOT run under Python 3 as cloned (checked 2026-09-21;
# there is no separate "Python 3 branch", contrary to some older advice). See the five patches
# in references/hess-local-h2.md's "Install" section -- apply them to your own clone once before
# using this script. HESS_DIR below must point at that patched clone.
#
# Usage:
#   bash hess_local_h2.sh <sumstats.txt> <chrom> <bfile_prefix> <partition.bed> <out_prefix>
#
# sumstats.txt columns (tab- or whitespace-delimited, header row): SNP CHR BP A1 A2 Z N
#   - CHR must match <chrom> exactly as a string (HESS does a literal string compare)
#   - A1/A2 must be single-letter alleles; SNP must be an rsID (HESS drops anything else)
# partition.bed columns (header row): chr  start  stop   (chr as "chr1", "chr2", ... "chrN")
# bfile_prefix: PLINK .bed/.bim/.fam reference panel covering the same chromosome/positions,
#   in-sample or ancestry-matched (HESS's quadratic form is unstable otherwise)

set -euo pipefail

SUMSTATS=${1:?Usage: hess_local_h2.sh <sumstats.txt> <chrom> <bfile_prefix> <partition.bed> <out_prefix>}
CHROM=${2:?provide chrom, e.g. 1}
BFILE=${3:?provide bfile_prefix}
PARTITION=${4:?provide partition.bed}
OUT=${5:?provide out_prefix}

HESS_DIR=${HESS_DIR:-./hess}

python "${HESS_DIR}/hess.py" \
    --local-hsqg "${SUMSTATS}" \
    --chrom "${CHROM}" \
    --bfile "${BFILE}" \
    --partition "${PARTITION}" \
    --out "${OUT}"

echo "Step 1 done: ${OUT}_chr${CHROM}.{info,eig,prjsq}.gz"
echo "Repeat for every chromosome you have data for, then run step 2 yourself once ALL 22"
echo "autosomes' step-1 output exists (HESS's own design requirement, not optional):"
echo "  python \${HESS_DIR}/hess.py --prefix ${OUT} --out <trait>_local_h2 --tot-hsqg <h2> <h2_SE>"
echo "(--tot-hsqg takes the genome-wide h2 estimate and SE from LDSC, used as the global constraint)"
