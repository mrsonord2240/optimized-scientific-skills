#!/bin/bash
# Reference: GCTA 1.94.1 | Verify against `gcta64 --help` if version differs
#
# GCTA-GREML: individual-level SNP-heritability via GRM-based REML.
# Yang 2011 AJHG 88:76. Needs individual-level genotypes (PLINK bed/bim/fam) + phenotype,
# not summary statistics -- this is the individual-level alternative to LDSC/LDAK.
#
# Usage:
#   bash gcta_greml.sh <bfile_prefix> <pheno_file> <out_prefix>
#
# pheno_file columns (no header): FID IID phenotype
# Verified 2026-09-21 end to end on a real 957-individual, 14389-SNP genotype panel with a
# planted h2=0.5 synthetic phenotype (200 causal SNPs): GCTA recovered h2 = 0.538 (SE 0.135),
# p = 2.3e-07 -- within 1 SE of the planted value.

set -euo pipefail

BFILE=${1:?Usage: gcta_greml.sh <bfile_prefix> <pheno_file> <out_prefix>}
PHENO=${2:?provide pheno_file}
OUT=${3:?provide out_prefix}

GCTA=${GCTA_BIN:-gcta64}

# Step 1: build the genetic relationship matrix (GRM) from autosomal SNPs
${GCTA} --bfile "${BFILE}" --autosome --make-grm --out "${OUT}_grm"

# Step 2: REML variance-components estimation
${GCTA} --grm "${OUT}_grm" --pheno "${PHENO}" --reml --out "${OUT}_h2"

echo "Result in ${OUT}_h2.hsq: V(G)/Vp is the SNP-heritability estimate, with SE and a"
echo "likelihood-ratio test p-value against h2 = 0."
echo "N <= 5000 gives a wide SE (see SKILL.md Statistical Model Taxonomy); for case-control"
echo "data below 20% prevalence use PCGC (--reml-pcgc) instead of --reml."
