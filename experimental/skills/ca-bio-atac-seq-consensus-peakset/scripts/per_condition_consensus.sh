#!/bin/bash
# Purpose: per-condition consensus (peak present in >= MIN_REPS replicates of that condition), then the union
#   across conditions, so condition-specific peaks survive.
# Usage:   per_condition_consensus.sh OUTDIR MIN_REPS cond1=rep1.narrowPeak,rep2.narrowPeak,rep3.narrowPeak cond2=...
# Output:  OUTDIR/<cond>_consensus.bed per condition and OUTDIR/all_conditions_consensus.bed (variable width).
# Next:    for differential counting re-center to fixed width, or use routes/iterative-overlap.md instead.
# Needs:   Bash, bedtools, awk, sort. Checked: not yet executed after the recut.
set -euo pipefail

if (($# < 3)); then
    echo "usage: $0 OUTDIR MIN_REPS cond=rep1,rep2[,rep3...] cond2=... (at least two conditions)" >&2
    exit 1
fi
OUTDIR=$1; MIN_REPS=$2; shift 2
mkdir -p "$OUTDIR"

for spec in "$@"; do
    cond=${spec%%=*}; files=${spec#*=}
    IFS=',' read -r -a reps <<< "$files"
    inputs=(); names=()
    for i in "${!reps[@]}"; do
        [[ -r "${reps[$i]}" ]] || { echo "ERROR: cannot read ${reps[$i]}" >&2; exit 1; }
        sort -k1,1 -k2,2n "${reps[$i]}" > "$OUTDIR/${cond}_rep$((i+1)).sorted.bed"
        inputs+=("$OUTDIR/${cond}_rep$((i+1)).sorted.bed"); names+=("rep$((i+1))")
    done
    bedtools multiinter -i "${inputs[@]}" -names "${names[@]}" | \
        awk -v m="$MIN_REPS" '$4 >= m {print $1"\t"$2"\t"$3}' > "$OUTDIR/${cond}_consensus.bed"
    echo "$cond: $(wc -l < "$OUTDIR/${cond}_consensus.bed") regions in >= $MIN_REPS of ${#reps[@]} replicates"
done

rm -f "$OUTDIR/all_conditions_consensus.bed"  # a rerun must not union its own earlier output
cat "$OUTDIR"/*_consensus.bed | sort -k1,1 -k2,2n | bedtools merge -i - > "$OUTDIR/all_conditions_consensus.bed"
echo "Union across conditions: $(wc -l < "$OUTDIR/all_conditions_consensus.bed") regions -> $OUTDIR/all_conditions_consensus.bed"
