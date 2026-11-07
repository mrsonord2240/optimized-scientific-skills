#!/bin/bash
# Checked on bcftools 1.19+, vcftools 0.1.16, peddy 0.4.8, somalier 0.2.19
# Sample identity QC: swaps, duplicates, relatedness, sex.
#
# Usage:
#   identity_qc.sh gtcheck  <in.vcf.gz> [truth.vcf.gz]        quick same-file swap/duplicate check; with a truth VCF, concordance to it
#   identity_qc.sh king     <in.vcf.gz>                       KING-robust kinship via vcftools --relatedness2
#   identity_qc.sh peddy    <in.vcf.gz> <cohort.ped>          sex, relatedness, ancestry vs the PED (human, genome-wide/exome only)
#   identity_qc.sh somalier <sites.vcf.gz> <ref.fa> <cohort.ped> <in.vcf.gz>   sketch extract then relate
# Env: OUT=<dir> (default identity_qc)   THREADS (peddy, default 4)

set -euo pipefail

MODE=${1:-}
OUT=${OUT:-identity_qc}
THREADS=${THREADS:-4}
mkdir -p "$OUT"

need() { for f in "$@"; do [ -f "$f" ] || { echo "ERROR: input not found: $f" >&2; exit 1; }; done; }

case "$MODE" in
    gtcheck)
        need "${2:-}"
        if [ -n "${3:-}" ]; then
            need "$3"
            bcftools gtcheck -g "$3" "$2" > "$OUT/gtcheck.txt"
        else
            bcftools gtcheck "$2" > "$OUT/gtcheck.txt"
        fi
        echo "DC lines = query, genotyped, discordance, sites: $OUT/gtcheck.txt" ;;
    king)
        need "${2:-}"
        vcftools --gzvcf "$2" --relatedness2 --out "$OUT/kin"
        echo "Done: $OUT/kin.relatedness2" ;;
    peddy)
        need "${2:-}" "${3:-}"
        python -m peddy -p "$THREADS" --plot --prefix "$OUT/cohort_qc" "$2" "$3"
        echo "Done: $OUT/cohort_qc.*" ;;
    somalier)
        need "${2:-}" "${3:-}" "${4:-}" "${5:-}"
        somalier extract -d "$OUT/extracted/" --sites "$2" -f "$3" "$5"
        somalier relate --ped "$4" -o "$OUT/somalier" "$OUT"/extracted/*.somalier
        echo "Done: $OUT/somalier.*" ;;
    *)
        echo "Usage: identity_qc.sh gtcheck|king|peddy|somalier ... (see header)" >&2
        exit 1 ;;
esac
