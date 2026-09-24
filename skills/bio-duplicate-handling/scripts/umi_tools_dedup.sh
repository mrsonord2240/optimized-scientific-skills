#!/bin/bash
# UMI-aware deduplication with umi_tools dedup (directional, seeded so the count is reproducible).
#
# Usage:  umi_tools_dedup.sh scrna       in.bam out.bam   # 10x / scRNA: cell barcode CB + UMI UB, coordinate-sorted + indexed input
#         umi_tools_dedup.sh bulk-paired in.bam out.bam   # bulk paired-end, UMI in the RX tag; sorts and indexes a copy first
# Needs:  umi_tools, samtools. Checked on umi_tools 1.1.6, samtools 1.24.
# Guards: scrna checks the CB tag exists (with absent CB/UB, --per-cell writes an EMPTY BAM and still exits 0);
#         both modes fail if the output has no records.

set -euo pipefail

MODE=${1:-}
IN=${2:-}
OUT=${3:-}
if [ -z "$MODE" ] || [ -z "$IN" ] || [ -z "$OUT" ]; then
    echo "Usage: umi_tools_dedup.sh scrna|bulk-paired in.bam out.bam" >&2
    exit 1
fi

case "$MODE" in
    scrna)
        # 10x / scRNA -- group by cell barcode + UMI. Check the tags exist first.
        n_cb=$(samtools view "$IN" | head -1000 | grep -c 'CB:Z:' || true)    # must be > 0
        if [ "$n_cb" -eq 0 ]; then
            echo "ERROR: no CB:Z: tags in the first 1000 records of $IN; scRNA dedup would write an empty BAM" >&2
            exit 2
        fi
        umi_tools dedup --stdin="$IN" --stdout="$OUT" \
            --extract-umi-method=tag --umi-tag=UB --cell-tag=CB \
            --per-cell --method=directional --random-seed=1
        ;;
    bulk-paired)
        # Bulk UMI, paired-end (UMI in the RX tag)
        TMP=$(mktemp -d)
        trap 'rm -rf "$TMP"' EXIT
        samtools sort -o "$TMP/sorted.bam" "$IN" && samtools index "$TMP/sorted.bam"
        umi_tools dedup --stdin="$TMP/sorted.bam" --stdout="$OUT" --paired \
            --extract-umi-method=tag --umi-tag=RX --method=directional --random-seed=1
        ;;
    *)
        echo "ERROR: mode must be scrna or bulk-paired" >&2
        exit 1 ;;
esac

n_out=$(samtools view -c "$OUT")
echo "dedup records: $n_out"
test "$n_out" -gt 0
