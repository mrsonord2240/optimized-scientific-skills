#!/usr/bin/env bash
# Plot the significant events of an rMATS event file with rmats2sashimiplot (checked on 4.0.0), one PDF
# per event per group, and fail if any figure is missing (rmats2sashimiplot exits 0 when it fails).
# Usage: rmats2sashimiplot_events.sh EVENTS.MATS.JC.txt EVENT_TYPE OUT_DIR B1_BAMS B2_BAMS [LABEL1 LABEL2]
#   EVENT_TYPE  SE | A5SS | A3SS | MXE | RI      B1_BAMS / B2_BAMS  comma-separated BAM lists (rMATS --b1 / --b2)
#   Env: FDR (0.05)  DPSI (0.1, |IncLevelDifference| cutoff)  COLORS ('#1f77b4,#ff7f0e', one per group)
# Writes OUT_DIR/ (the rmats2sashimiplot output), OUT_DIR_sig.txt (filtered events) and OUT_DIR_grouping.gf.
set -euo pipefail
[ $# -ge 5 ] || { sed -n '2,7p' "$0" >&2; exit 2; }
events=$1; etype=$2; out=${3%/}; b1=$4; b2=$5; l1=${6:-Control}; l2=${7:-Treatment}
FDR=${FDR:-0.05}; DPSI=${DPSI:-0.1}; COLORS=${COLORS:-#1f77b4,#ff7f0e}
n1=$(echo "$b1" | tr ',' '\n' | wc -l); n2=$(echo "$b2" | tr ',' '\n' | wc -l)

# rmats2sashimiplot plots every row: keep only significant events (columns found by header name)
awk -F'\t' -v fdr="$FDR" -v dpsi="$DPSI" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; print; next}
    $c["FDR"]<fdr && ($c["IncLevelDifference"]>dpsi || $c["IncLevelDifference"]<-dpsi)' \
    "$events" > "${out}_sig.txt"
n_events=$(( $(wc -l < "${out}_sig.txt") - 1 ))
[ "$n_events" -gt 0 ] || { echo "no event passes FDR<$FDR and |dPSI|>$DPSI in $events; nothing plotted" >&2; exit 0; }

# group file: "label: first-last", 1-based over the --b1 replicates then the --b2 replicates
printf '%s: 1-%d\n%s: %d-%d\n' "$l1" "$n1" "$l2" $((n1 + 1)) $((n1 + n2)) > "${out}_grouping.gf"

rmats2sashimiplot \
    --b1 "$b1" --b2 "$b2" \
    --event-type "$etype" \
    -e "${out}_sig.txt" \
    --l1 "$l1" --l2 "$l2" \
    -o "$out" \
    --exon_s 1 --intron_s 5 \
    --group-info "${out}_grouping.gf" \
    --color "$COLORS"

n_pdf=$(find "$out/Sashimi_plot" -name '*.pdf' -size +0 2>/dev/null | wc -l)
[ "$n_pdf" -eq "$n_events" ] || { echo "rmats2sashimiplot wrote $n_pdf of $n_events figures" >&2; exit 1; }
echo "OK: $n_pdf figures in $out/Sashimi_plot"
