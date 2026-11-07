#!/bin/bash
# Checked on bcftools 1.24
# Normalize two or more VCFs the same way, then intersect them with bcftools isec.
#
# Usage: [ATOMIZE=0] compare_callers.sh <reference.fa> <outdir> <a.vcf.gz> <b.vcf.gz> [...]
# isec output in <outdir>/isec for two inputs: 0000.vcf private to the first, 0001.vcf private to the second,
# 0002.vcf / 0003.vcf shared (records from each file).

set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REF=${1:-}
OUT=${2:-}
if [ -z "$REF" ] || [ -z "$OUT" ] || [ "$#" -lt 4 ]; then
    echo "Usage: compare_callers.sh reference.fa outdir a.vcf.gz b.vcf.gz [...]" >&2
    exit 1
fi
shift 2
mkdir -p "$OUT"

NORMED=()
for vcf in "$@"; do
    base=$(basename "$vcf" .vcf.gz)
    bash "$HERE/normalize_vcf.sh" "$REF" "$vcf" "$OUT/$base.norm.vcf.gz" > "$OUT/$base.norm.log"
    tail -n 8 "$OUT/$base.norm.log" | head -n 5
    NORMED+=("$OUT/$base.norm.vcf.gz")
done

bcftools isec -p "$OUT/isec" "${NORMED[@]}"
echo "Done: $OUT/isec"
