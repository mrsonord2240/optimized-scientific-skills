#!/bin/bash
# Reference: pysam 0.24.1, samtools 1.24 (checked); pysam 0.22+, samtools 1.19+
# Convert between SAM/BAM/CRAM formats

set -e

INPUT=$1
OUTPUT=$2
REFERENCE=$3

if [ -z "$INPUT" ] || [ -z "$OUTPUT" ]; then
    echo "Usage: convert_formats.sh <input> <output> [reference.fa]"
    echo "Examples:"
    echo "  convert_formats.sh input.sam output.bam"
    echo "  convert_formats.sh input.bam output.cram reference.fa"
    echo "  convert_formats.sh input.cram output.bam reference.fa"
    exit 1
fi

if [ "$INPUT" -ef "$OUTPUT" ]; then
    echo "Error: input and output are the same file; refusing to overwrite the input"
    exit 1
fi

EXT=$(echo "${OUTPUT##*.}" | tr '[:upper:]' '[:lower:]')

case "$EXT" in
    sam|bam) ;;
    cram)
        if [ -z "$REFERENCE" ]; then
            echo "Error: CRAM conversion requires reference.fa"
            exit 1
        fi
        ;;
    *)
        echo "Unknown output format: $EXT"
        exit 1
        ;;
esac

# -T is needed for CRAM output and for CRAM input, and is ignored for SAM/BAM
REF_ARGS=()
if [ -n "$REFERENCE" ]; then
    REF_ARGS=(-T "$REFERENCE")
fi

# From here on samtools writes OUTPUT: a failed run must not leave a header-only or partial file behind
trap 'rc=$?; if [ "$rc" -ne 0 ]; then rm -f "$OUTPUT"; fi' EXIT

case "$EXT" in
    sam) samtools view -h "${REF_ARGS[@]}" -o "$OUTPUT" "$INPUT" ;;
    bam) samtools view -b "${REF_ARGS[@]}" -o "$OUTPUT" "$INPUT" ;;
    cram) samtools view -C "${REF_ARGS[@]}" -o "$OUTPUT" "$INPUT" ;;
esac

echo "Converted $INPUT -> $OUTPUT"
