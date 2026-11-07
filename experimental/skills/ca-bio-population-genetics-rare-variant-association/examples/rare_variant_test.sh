#!/usr/bin/env bash
# Reference: regenie 4.1.3 | Verify API if version differs
# Gene/region-based rare-variant aggregation with regenie: build masks (functional class x MAF
# cutoff) and run burden plus SKAT-O and the ACAT-O omnibus, with Firth keeping an imbalanced
# binary trait calibrated. Reuses a step-1 whole-genome ridge null (LOCO) so the relatedness model
# matches single-variant GWAS. The mask IS the hypothesis: this runs LoF and LoF+missense masks at
# two nested MAF cutoffs and lets ACAT-O combine them. Outputs go to a caller-supplied directory
# (default: a fresh temp dir) so nothing is written to the current dir.
# Step 1 and step 2 take DIFFERENT genotype sets: step 1 needs QC'd common (array-like) variants
# (rare variants stop the ridge fit with "low variance"); step 2 takes the exome/rare-variant set.
# Usage: ./rare_variant_test.sh <array_prefix> <exome_prefix> <pheno.txt> <covar.txt> <anno.txt> <sets.txt> <masks.txt> [output_dir]
set -euo pipefail

USAGE="Usage: $0 <array_prefix> <exome_prefix> <pheno.txt> <covar.txt> <anno.txt> <sets.txt> <masks.txt> [output_dir]"
ARRAY="${1:?$USAGE}"
EXOME="${2:?$USAGE}"
PHENO="${3:?pheno file required}"
COVAR="${4:?covar file required}"
ANNO="${5:?annotation file required}"
SETS="${6:?set-list file required}"
MASKS="${7:?mask-definition file required}"
OUTDIR="${8:-$(mktemp -d)}"
mkdir -p "$OUTDIR"
echo "Outputs -> $OUTDIR"

# Validate the mask inputs FIRST. --check-burden-files reports set-list variants missing from the
# annotation file - a silent source of empty or wrong masks - into <out>_masks_report.txt.
# --ignore-pred is required: this validation runs before step 1, so no LOCO predictor exists yet.
regenie --step 2 --bed "$EXOME" --phenoFile "$PHENO" --covarFile "$COVAR" \
    --anno-file "$ANNO" --set-list "$SETS" --mask-def "$MASKS" \
    --aaf-bins 0.001,0.01 --build-mask max --check-burden-files --ignore-pred \
    --bt --out "$OUTDIR/check"

# Step 1 builds the LOCO whole-genome ridge predictor (the null) once, on the common-variant set.
# --bsize 1000 is the SNP block size for the ridge stacking; --lowmem spills to disk to bound memory.
regenie --step 1 --bed "$ARRAY" --phenoFile "$PHENO" --covarFile "$COVAR" \
    --bsize 1000 --bt --lowmem --lowmem-prefix "$OUTDIR/tmp_rg" --out "$OUTDIR/fit_null"

# Step 2 set tests on the exome set. --aaf-bins sets the MAF ceilings (an all-variant ".all" mask is
# added on top; singleton masks need --singleton-carrier).
# --vc-tests skato,acato requests the SKAT-O test and the ACAT omnibus alongside the default burden.
# --firth --approx gives the penalized-LRT fallback that keeps an imbalanced binary tail calibrated;
# --pThresh 0.05 restricts the (slower) Firth step to masks below nominal significance.
regenie --step 2 --bed "$EXOME" --phenoFile "$PHENO" --covarFile "$COVAR" \
    --pred "$OUTDIR/fit_null_pred.list" \
    --anno-file "$ANNO" --set-list "$SETS" --mask-def "$MASKS" \
    --aaf-bins 0.001,0.01 --build-mask max --vc-tests skato,acato \
    --bt --firth --approx --pThresh 0.05 --out "$OUTDIR/gene_tests"

echo "Gene-test results: $OUTDIR/gene_tests_*.regenie"
echo "Mask concordance report: $OUTDIR/check_masks_report.txt"
# Exome-wide gene significance ~2.5e-6 (Bonferroni 0.05 over ~20,000 genes); tighten further across
# the multiple masks tested here, or rely on the ACAT-O omnibus to absorb them.
