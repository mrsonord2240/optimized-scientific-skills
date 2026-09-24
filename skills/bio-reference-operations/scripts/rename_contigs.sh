#!/bin/bash
# Rename BAM contigs (@SQ SN) without re-aligning: header-only change via `samtools reheader`.
# Only rename when the sequences are the same (compare LN, and M5 where the header has it).
# Inputs:  in.bam; map.tsv = old<TAB>new per line (e.g. chr22<TAB>22). Contigs not in the map keep their names.
#          Build map.tsv from an NCBI assembly report with the awk line in references/contig-naming.md.
# Output:  out.bam (+ .bai). A failed reheader leaves no out.bam. Ends with a name/length check against ref.fa.fai when given.
# Usage:   rename_contigs.sh in.bam map.tsv out.bam [ref.fa]     (ref.fa.fai must exist for the check)
# Tested:  samtools 1.24
set -euo pipefail

IN=${1:?usage: rename_contigs.sh in.bam map.tsv out.bam [ref.fa]}
MAP=${2:?map.tsv}
OUT=${3:?out.bam}
REF=${4:-}

HDR="${OUT}.hdr"
samtools view -H "$IN" | awk -F'\t' -v OFS='\t' 'NR==FNR{m[$1]=$2; next}
    /^@SQ/{for(i=2;i<=NF;i++) if($i~/^SN:/){n=substr($i,4); if(n in m) $i="SN:" m[n]}} {print}' "$MAP" - > "$HDR"
samtools reheader "$HDR" "$IN" > "${OUT}.tmp" && mv "${OUT}.tmp" "$OUT" && samtools index "$OUT"
rm -f "$HDR"

if [ -n "$REF" ]; then
    diff <(samtools view -H "$OUT" | awk '/^@SQ/{print $2, $3}') <(awk '{print "SN:"$1, "LN:"$2}' "${REF}.fai") && echo OK
fi
