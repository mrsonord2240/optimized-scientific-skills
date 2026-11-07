#!/bin/bash
# Checked on bcftools 1.24
# Count REF alleles that do not match the reference, without changing anything.
#
# Usage: check_ref.sh <reference.fa> <input.vcf.gz>
# Exit 0 if none, 1 if any. The full list: bcftools norm -f ref.fa -c w input.vcf.gz -Ou > /dev/null

set -euo pipefail

REF=${1:-}
VCF=${2:-}
if [ -z "$REF" ] || [ -z "$VCF" ]; then
    echo "Usage: check_ref.sh reference.fa input.vcf.gz" >&2
    exit 1
fi
for f in "$REF" "$VCF"; do
    [ -f "$f" ] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

N=$(bcftools norm -f "$REF" -c w "$VCF" -Ou 2>&1 >/dev/null | grep -c 'REF_MISMATCH\|does not match' || true)
if [ "$N" -gt 0 ]; then
    echo "$N REF allele(s) do not match $REF (wrong build, contig naming, or masking?)"
    exit 1
fi
echo "REF alleles match $REF"
