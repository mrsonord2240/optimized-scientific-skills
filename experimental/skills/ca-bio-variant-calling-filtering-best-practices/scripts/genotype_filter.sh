#!/bin/bash
# Genotype-level filter: set untrustworthy per-sample genotypes to no-call (./.) and keep the site.
# Run AFTER the site filter and BEFORE any cohort missingness, HWE or allele-frequency work.
#
# Usage:  genotype_filter.sh <in.vcf[.gz]> <out.vcf.gz> [gq=20] [dp=8] [ab=0]
#           gq=20   set ./. where FMT/GQ < gq (0 disables)
#           dp=8    set ./. where FMT/DP < dp (0 disables; 8-10 for WGS)
#           ab=1    also set ./. at hets whose alt fraction (from AD) is <= 0.2 or >= 0.8 (default off)
# Needs:  bcftools. Checked on bcftools 1.24.
# Traps:  exit 2 when FORMAT/GQ, DP or AD is absent from the header for an enabled filter; always prints how many
#         genotypes became no-call, so a setting that changes nothing is visible.

set -euo pipefail

IN=${1:-}
OUT=${2:-}
if [ -z "$IN" ] || [ -z "$OUT" ]; then
    echo "Usage: genotype_filter.sh <in.vcf[.gz]> <out.vcf.gz> [gq=20] [dp=8] [ab=0]" >&2
    exit 1
fi
shift 2
GQ=20; DP=8; AB=0
for kv in "$@"; do
    case "$kv" in
        gq=*) GQ=${kv#gq=} ;;
        dp=*) DP=${kv#dp=} ;;
        ab=*) AB=${kv#ab=} ;;
        *) echo "ERROR: unknown argument $kv" >&2; exit 1 ;;
    esac
done

HDR=$(bcftools view -h "$IN")
need() {
    if ! grep -q "^##FORMAT=<ID=$1," <<<"$HDR"; then
        echo "ERROR: FORMAT/$1 is not in the header of $IN; cannot apply that filter (set $2=0 to skip it)." >&2
        exit 2
    fi
}
EXPR=
add() { if [ -z "$EXPR" ]; then EXPR="$1"; else EXPR="$EXPR | $1"; fi; }
if [ "$GQ" != 0 ]; then need GQ gq; add "FMT/GQ<$GQ"; fi
if [ "$DP" != 0 ]; then need DP dp; add "FMT/DP<$DP"; fi
if [ "$AB" = 1 ]; then
    need AD ab
    add 'GT="het" & (FMT/AD[:1]/(FMT/AD[:0]+FMT/AD[:1])<=0.2 | FMT/AD[:1]/(FMT/AD[:0]+FMT/AD[:1])>=0.8)'
fi

count_missing() { bcftools query -f '[%GT\n]' "$1" | awk '$0 ~ /\./ {n++} END{print n+0}'; }
before=$(count_missing "$IN")
if [ -z "$EXPR" ]; then
    echo "No filter enabled (gq=0 dp=0 ab=0): copying input unchanged."
    bcftools view "$IN" -Oz -o "$OUT"
else
    bcftools filter -S . -e "$EXPR" "$IN" -Oz -o "$OUT"
fi
bcftools index -t -f "$OUT"
after=$(count_missing "$OUT")
echo "expression: ${EXPR:-none}"
echo "no-call genotypes: before $before, after $after, newly set $((after - before))"
echo "records: $(bcftools view -H "$OUT" | wc -l)"
