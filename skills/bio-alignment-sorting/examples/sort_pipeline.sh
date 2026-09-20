#!/bin/bash
# Reference: bwa 0.7.17+, samtools 1.19+ (checked on bwa 0.7.19, samtools 1.24) | Verify API if version differs
# Align, sort, and index reads. Prerequisite: `bwa index <reference.fa>` has been run once.

# pipefail: without it a failing `bwa mem` in the pipe still gives exit 0 and a truncated BAM
set -eo pipefail

REF=$1
R1=$2
OUTPUT=$3
R2=$4
THREADS=${5:-8}

if [ -z "$REF" ] || [ -z "$R1" ] || [ -z "$OUTPUT" ]; then
    echo "Usage: sort_pipeline.sh <reference.fa> <R1.fq> <output.bam> [R2.fq] [threads]"
    echo "Example (paired): sort_pipeline.sh ref.fa reads_R1.fq aligned.bam reads_R2.fq 8"
    echo "Example (single): sort_pipeline.sh ref.fa reads.fq aligned.bam"
    exit 1
fi

for f in "$REF" "$R1" ${R2:+"$R2"}; do
    [ -f "$f" ] || { echo "ERROR: input not found: $f" >&2; exit 1; }
done
[ -f "$REF.bwt" ] || { echo "ERROR: no bwa index for $REF; run: bwa index $REF" >&2; exit 1; }

echo "Aligning and sorting..."

if [ -z "$R2" ]; then
    # Single-end (R2 omitted)
    bwa mem -t "$THREADS" "$REF" "$R1" | samtools sort -@ 4 -o "$OUTPUT"
    EXPECTED=$(( $(gzip -cdf "$R1" | wc -l) / 4 ))
else
    # Paired-end
    bwa mem -t "$THREADS" "$REF" "$R1" "$R2" | samtools sort -@ 4 -o "$OUTPUT"
    EXPECTED=$(( ($(gzip -cdf "$R1" | wc -l) + $(gzip -cdf "$R2" | wc -l)) / 4 ))
fi

# bwa can exit 0 on a truncated FASTQ; every input read must appear once as a primary record
GOT=$(samtools view -c -F 0x900 "$OUTPUT")
if [ "$GOT" -ne "$EXPECTED" ]; then
    echo "ERROR: $OUTPUT has $GOT primary records, expected $EXPECTED input reads" >&2
    exit 1
fi

echo "Indexing..."
samtools index "$OUTPUT"

echo "Done: $OUTPUT"
samtools flagstat "$OUTPUT"
