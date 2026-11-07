#!/bin/bash
# Checked on samtools 1.24
# Sort a BAM/CRAM and check the result.
#
# Usage: sort_bam.sh <order> <input> <output> [threads] [mem_per_thread]
#   order: coordinate | name-natural | name-ascii | tag:<TAG> | template-coordinate
# Env:   REF       reference FASTA, required when <output> ends in .cram (and for a CRAM input)
#        SORT_TMP  temp-file prefix on fast local storage (default: a fresh temp directory)

set -euo pipefail

ORDER=${1:-}
INPUT=${2:-}
OUTPUT=${3:-}
THREADS=${4:-4}
MEM=${5:-2G}

if [ -z "$ORDER" ] || [ -z "$INPUT" ] || [ -z "$OUTPUT" ]; then
    echo "Usage: sort_bam.sh <coordinate|name-natural|name-ascii|tag:TAG|template-coordinate> in.bam out.bam [threads] [mem_per_thread]" >&2
    echo "       name-natural: only samtools reads the result; name-ascii: Picard, GATK or htsjdk reads it (routes/required-order.md)" >&2
    exit 1
fi
[ -f "$INPUT" ] || { echo "ERROR: input not found: $INPUT" >&2; exit 1; }

case "$ORDER" in
    coordinate)           FLAGS=() ;;
    name-natural)         FLAGS=(-n) ;;
    name-ascii)           FLAGS=(-N) ;;
    tag:?*)               FLAGS=(-t "${ORDER#tag:}") ;;
    template-coordinate)  FLAGS=(--template-coordinate) ;;
    *) echo "ERROR: unknown order '$ORDER'. Use coordinate, name-natural, name-ascii, tag:<TAG> or template-coordinate." >&2; exit 2 ;;
esac

OUTFMT=bam
REFARGS=()
case "$OUTPUT" in *.cram) OUTFMT=cram ;; esac
if [ "$OUTFMT" = cram ] || [ -n "${REF:-}" ]; then
    [ -n "${REF:-}" ] || { echo "ERROR: CRAM output needs REF=<reference.fa>" >&2; exit 2; }
    REFARGS=(--reference "$REF")
fi

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
PREFIX=${SORT_TMP:-$TMP/sort}
PART="$TMP/out.$OUTFMT"

echo "Sorting $INPUT ($ORDER, $THREADS threads, -m $MEM)..."
if [ "$ORDER" = template-coordinate ]; then
    # The sort needs MC tags; fixmate -m adds them and needs mates adjacent.
    samtools collate -O -u -@ "$THREADS" "${REFARGS[@]}" "$INPUT" "$TMP/collate" | \
        samtools fixmate -m -u -@ "$THREADS" - - | \
        samtools sort "${FLAGS[@]}" -@ "$THREADS" -m "$MEM" -T "$PREFIX" -O "$OUTFMT" "${REFARGS[@]}" -o "$PART" -
else
    samtools sort "${FLAGS[@]}" -@ "$THREADS" -m "$MEM" -T "$PREFIX" -O "$OUTFMT" "${REFARGS[@]}" -o "$PART" "$INPUT"
fi

# rc 16 = truncated, rc 4 = zero-byte
samtools quickcheck -v "$PART" || { echo "ERROR: sorted output failed quickcheck; input kept, output not published" >&2; exit 3; }

mv "$PART" "$OUTPUT"
echo "Header: $(samtools view -H "$OUTPUT" | grep '^@HD' || echo 'no @HD line')"

if [ "$ORDER" = coordinate ]; then
    echo "Indexing..."
    samtools index "$OUTPUT"
    echo "Indexed: $OUTPUT"
else
    echo "Not indexable (only coordinate-sorted files can be indexed)."
fi
echo "Done: $OUTPUT"
