#!/bin/bash
# Downsample a BAM (hash-based `samtools view -s`, pair-consistent) to a target read count, or to
# the read count of another BAM (tumor-normal matching). Lands a few % off the target.
# Never upsamples: a fraction of 1 or more spliced into `-s` (`1.084419`) silently keeps only 8%,
# so a BAM already at or under the target is copied unchanged.
#
# Usage: match_read_count.sh input.bam output.bam --target N       # N primary reads (counts -F 2304)
#        match_read_count.sh input.bam output.bam --like other.bam # other.bam's mapped primary reads (-F 2308)
# Checked on samtools 1.24.
set -euo pipefail

input=${1:?usage: $0 input.bam output.bam --target N | --like other.bam}
output=${2:?output.bam}
mode=${3:?--target N | --like other.bam}
value=${4:?value}

case "$mode" in
    --target)
        total=$(samtools view -c -F 2304 "$input")
        target=$value
        ;;
    --like)
        total=$(samtools view -c -F 2308 "$input")
        target=$(samtools view -c -F 2308 "$value")
        ;;
    *) echo "unknown mode $mode" >&2; exit 2 ;;
esac

if [ "$total" -le "$target" ]; then
    echo "only $total reads, not more than the target $target; copying unchanged" >&2
    cp "$input" "$output"
else
    frac=$(awk -v t="$target" -v n="$total" 'BEGIN{printf "%.6f", t/n}')
    samtools view -s "1.${frac#*.}" -b -o "$output" "$input"
fi
