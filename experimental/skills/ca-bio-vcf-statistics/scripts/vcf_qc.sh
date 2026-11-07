#!/bin/bash
# Checked on bcftools 1.19+
# Callset-level and per-sample VCF statistics with bcftools stats.
#
# Usage: vcf_qc.sh <input.vcf.gz|bcf> <outdir>
# Env:   PLOT=1   also run plot-vcfstats (needs matplotlib; summary.pdf needs pdflatex or tectonic)

set -euo pipefail

VCF=${1:-}
OUT=${2:-}
if [ -z "$VCF" ] || [ -z "$OUT" ]; then
    echo "Usage: [PLOT=1] vcf_qc.sh input.vcf.gz outdir" >&2
    exit 1
fi
[ -f "$VCF" ] || { echo "ERROR: input not found: $VCF" >&2; exit 1; }
mkdir -p "$OUT"

bcftools stats "$VCF" > "$OUT/stats.txt"
bcftools stats -s - "$VCF" > "$OUT/per_sample.txt"

{
    echo "== SN: summary numbers"
    grep '^SN' "$OUT/stats.txt" | cut -f3-
    echo
    echo "== TSTV: ts, tv, ts/tv (sites; the last field is ts/tv of singletons)"
    grep '^TSTV' "$OUT/stats.txt" | cut -f2-
    echo
    echo "== PSC: per-sample [3]sample [4]nRefHom [5]nNonRefHom [6]nHets [7]nTransitions [8]nTransversions [9]nIndels [10]avgDP [11]nSingletons [12]nMissing"
    grep '^PSC' "$OUT/per_sample.txt" | cut -f3-
} > "$OUT/summary.txt"

N_ALL=$(bcftools view -H "$VCF" | wc -l)
N_SNP=$(bcftools view -v snps -H "$VCF" | wc -l)
N_INDEL=$(bcftools view -v indels -H "$VCF" | wc -l)
N_PASS=$(bcftools view -f PASS -H "$VCF" | wc -l)
N_NOTFAIL=$(bcftools view -f .,PASS -H "$VCF" | wc -l)
MEANQ=$(bcftools query -f '%QUAL\n' "$VCF" | awk '$1!="." {s+=$1; n++} END{if(n) printf "%.1f", s/n; else print "NA"}')
{
    echo "records=$N_ALL snps=$N_SNP indels=$N_INDEL"
    echo "FILTER=PASS strictly: $N_PASS; not failed (PASS or '.'): $N_NOTFAIL (use the second on an unfiltered callset)"
    echo "mean QUAL=$MEANQ"
} | tee "$OUT/counts.txt"

cat "$OUT/summary.txt"

if [ "${PLOT:-0}" = 1 ]; then
    plot-vcfstats -p "$OUT/plots/" "$OUT/stats.txt" || echo "WARNING: plot-vcfstats failed (matplotlib, or no pdflatex/tectonic for summary.pdf; PNGs may still exist)" >&2
fi
echo "Done: $OUT/stats.txt per_sample.txt summary.txt counts.txt"
