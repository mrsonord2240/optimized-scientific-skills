#!/bin/bash
# Checked on bcftools 1.19+
# bcftools stats per region class (e.g. GIAB stratification BEDs), with Ti/Tv and SNP/indel counts side by side.
#
# Usage: stratified_stats.sh <input.vcf.gz> <outdir> <name=regions.bed> [<name=regions.bed> ...]
# The input must be bgzipped and indexed (-R needs the index).

set -euo pipefail

VCF=${1:-}
OUT=${2:-}
if [ -z "$VCF" ] || [ -z "$OUT" ] || [ "$#" -lt 3 ]; then
    echo "Usage: stratified_stats.sh input.vcf.gz outdir easy=easy.bed difficult=difficult.bed ..." >&2
    exit 1
fi
shift 2
[ -f "$VCF" ] || { echo "ERROR: input not found: $VCF" >&2; exit 1; }
mkdir -p "$OUT"

printf 'stratum\tsnps\tindels\tTi/Tv\n'
for spec in "$@"; do
    name=${spec%%=*}
    bed=${spec#*=}
    [ -f "$bed" ] || { echo "ERROR: BED not found: $bed" >&2; exit 1; }
    bcftools stats -R "$bed" "$VCF" > "$OUT/$name.txt"
    snps=$(grep '^SN' "$OUT/$name.txt" | grep 'number of SNPs:' | cut -f4)
    indels=$(grep '^SN' "$OUT/$name.txt" | grep 'number of indels:' | cut -f4)
    tstv=$(grep '^TSTV' "$OUT/$name.txt" | cut -f5)
    printf '%s\t%s\t%s\t%s\n' "$name" "$snps" "$indels" "$tstv"
done
