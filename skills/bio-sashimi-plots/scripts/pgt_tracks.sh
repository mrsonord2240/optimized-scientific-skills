#!/usr/bin/env bash
# Inputs for a pyGenomeTracks figure: one coverage bedGraph per group and a BEDPE of junction arcs (checked on
# samtools 1.24, bedtools 2.31, regtools 1.0.0).
# Usage: pgt_tracks.sh OUT_DIR CTRL_BAMS TRT_BAMS      (BAM lists comma-separated, each BAM indexed)
# Writes OUT_DIR/{ctrl,trt}.bedgraph, OUT_DIR/junctions.bedpe (score = reads summed over all BAMs; `-s XS` needs XS-tagged BAMs,
# otherwise the strand is `?`) and the merged BAMs.
set -euo pipefail
[ $# -eq 3 ] || { sed -n '2,7p' "$0" >&2; exit 2; }
out=${1%/}; IFS=, read -r -a ctrl <<< "$2"; IFS=, read -r -a trt <<< "$3"
mkdir -p "$out"

# one merged BAM per group; -split is essential: without it introns are filled with coverage
samtools merge -f "$out/ctrl_merged.bam" "${ctrl[@]}" && samtools index "$out/ctrl_merged.bam"
samtools merge -f "$out/trt_merged.bam" "${trt[@]}" && samtools index "$out/trt_merged.bam"
bedtools genomecov -ibam "$out/ctrl_merged.bam" -split -bga > "$out/ctrl.bedgraph"
bedtools genomecov -ibam "$out/trt_merged.bam" -split -bga > "$out/trt.bedgraph"

samtools merge -f "$out/all_merged.bam" "$out/ctrl_merged.bam" "$out/trt_merged.bam" && samtools index "$out/all_merged.bam"   # regtools needs an indexed BAM
regtools junctions extract -s XS -o "$out/regtools_junctions.bed" "$out/all_merged.bam"
# regtools BED12 column 11 is blockSizes (anchor_left, anchor_right);
# column 12 is blockStarts (0, intron_length + anchor_left).
# Intron start = chromStart + anchor_left = $2 + a[1]
# Intron end   = chromStart + blockStarts[2] = $2 + b[2]
awk 'BEGIN{OFS="\t"} {split($11,a,","); split($12,b,","); s=$2+a[1]; e=$2+b[2]; print $1, s, s+1, $1, e-1, e, $5}' \
    "$out/regtools_junctions.bed" > "$out/junctions.bedpe"
[ -s "$out/junctions.bedpe" ] || { echo "no junctions written (unindexed BAM, or no split reads)" >&2; exit 1; }
echo "OK: $(wc -l < "$out/junctions.bedpe") junctions"
