#!/usr/bin/env bash
# Jutils on rMATS output: convert, heatmap, sashimi and a JC-vs-JCEC Venn (checked on Jutils 1.5).
# Usage: jutils_pipeline.sh JUTILS_DIR RMATS_DIR META.tsv BAM_LIST.tsv ANNOTATION.gtf COORDINATE OUT_DIR
#   META.tsv sample<TAB>condition   BAM_LIST.tsv sample<TAB>bam<TAB>condition   COORDINATE chr:start-end
#   The heatmap needs >= 2 events passing --q-value (Env Q, default 0.05); with fewer it prints "Skipping" and writes nothing.
# Writes OUT_DIR/{jutils_out,hm,sh,vn}. This script covers rMATS input; the leafcutter/MntJULiP/MAJIQ
# converters (--leafcutter-dir/--mntjulip-dir/--majiq-dir) are documented and verified in references/jutils.md.
set -euo pipefail
[ $# -eq 7 ] || { sed -n '2,6p' "$0" >&2; exit 2; }
ju=$(cd "$1" && pwd)/jutils.py; rmats=$2; meta=$3; bams=$4; gtf=$5; coord=$6; out=$7; Q=${Q:-0.05}
mkdir -p "$out"

# writes rmats_JC_results.tsv and rmats_JCEC_results.tsv into --out-dir
python3 "$ju" convert-results --rmats-dir "$rmats" --out-dir "$out/jutils_out/"

# Needs >= 2 events passing the cutoffs; writes clustermap*.pdf
python3 "$ju" heatmap --tsv-file "$out/jutils_out/rmats_JC_results.tsv" --meta-file "$meta" --q-value "$Q" --out-dir "$out/hm/" --pdf

python3 "$ju" sashimi --tsv-file "$out/jutils_out/rmats_JC_results.tsv" --meta-file "$meta" \
    --gtf "$gtf" --coordinate "$coord" --bam-list "$bams" --out-dir "$out/sh/" --pdf

# --tsv-file-list is a FILE with one "path<TAB>label" line per TSV, not a comma-separated list
printf '%s\trMATS_JC\n%s\trMATS_JCEC\n' "$out/jutils_out/rmats_JC_results.tsv" "$out/jutils_out/rmats_JCEC_results.tsv" > "$out/tsv_list.txt"
python3 "$ju" venn-diagram --tsv-file-list "$out/tsv_list.txt" --out-dir "$out/vn/"
