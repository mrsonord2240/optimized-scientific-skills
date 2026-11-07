#!/bin/bash
# Checked on vcftools 0.1.16
# Per-sample and per-site missingness, and exact-test HWE per site.
#
# Usage: missing_hwe.sh <input.vcf.gz> <outdir> [controls.vcf.gz]
#   HWE runs on controls.vcf.gz when given, else on the input.
# Apply genotype-level filters (set ./. where GQ/DP fail) BEFORE this; run it on one ancestry-homogeneous
# group at a time. Sample and variant missingness are coupled: drop the worst samples, re-run, then the worst variants.

set -euo pipefail

VCF=${1:-}
OUT=${2:-}
CTRL=${3:-$VCF}
if [ -z "$VCF" ] || [ -z "$OUT" ]; then
    echo "Usage: missing_hwe.sh input.vcf.gz outdir [controls.vcf.gz]" >&2
    exit 1
fi
for f in "$VCF" "$CTRL"; do
    [ -f "$f" ] || { echo "ERROR: input not found: $f" >&2; exit 1; }
done
mkdir -p "$OUT"

vcftools --gzvcf "$VCF" --missing-indv --out "$OUT/sample_miss"
vcftools --gzvcf "$VCF" --missing-site --out "$OUT/site_miss"
vcftools --gzvcf "$CTRL" --hardy --out "$OUT/hwe"

echo "Worst 10 samples by F_MISS:"
tail -n +2 "$OUT/sample_miss.imiss" | sort -k5,5gr | head -10
echo "Done: $OUT/sample_miss.imiss site_miss.lmiss hwe.hwe"
