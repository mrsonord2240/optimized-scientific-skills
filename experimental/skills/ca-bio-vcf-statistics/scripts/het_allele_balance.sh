#!/bin/bash
# Checked on bcftools 1.19+
# Mean allele balance of biallelic 0/1 hets per sample, from FORMAT/AD. Expect ~0.5 per sample.
# A shift off 0.5 with a high het count is a contamination signal, not the measurement.
#
# Usage: het_allele_balance.sh <input.vcf.gz>
# Handles 0/1, 1/0, 0|1, 1|0 only; use a multiallelic-aware parser for other genotypes.

set -euo pipefail

VCF=${1:-}
[ -n "$VCF" ] || { echo "Usage: het_allele_balance.sh input.vcf.gz" >&2; exit 1; }
[ -f "$VCF" ] || { echo "ERROR: input not found: $VCF" >&2; exit 1; }

bcftools query -f '[%SAMPLE\t%GT\t%AD\n]' "$VCF" | \
    awk -F'\t' '($2 ~ /^(0[\/|]1|1[\/|]0)$/) {split($3,a,","); d=a[1]+a[2]; if (d>0) {s[$1]+=a[2]/d; n[$1]++}}
        END {for (k in s) printf "%s\tmean het AB: %.3f (n=%d)\n", k, s[k]/n[k], n[k]}'
