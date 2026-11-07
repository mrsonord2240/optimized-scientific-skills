#!/bin/bash
# Reference: GATK 4.5+, bcftools 1.19+, pysam 0.22+, samtools 1.19+ | Verify API if version differs
# Checked on samtools 1.24, GATK 4.6.2.0, Picard 3.5.0
# Prepare reference genome for analysis

set -e

REF=$1

if [ -z "$REF" ]; then
    echo "Usage: prepare_reference.sh <reference.fa|.fasta|.fna[.gz]>  (.gz must be bgzip, not plain gzip)"
    exit 1
fi

if [ ! -f "$REF" ]; then
    echo "Error: Reference file not found: $REF"
    exit 1
fi

# GATK requires <name>.dict: drop .gz, then the last extension (genome.fasta -> genome.dict)
# GATK ignores genome.fasta.dict ("Fasta dict file genome.dict ... does not exist"); Picard accepts both names
NAME=$(basename "$REF")
NAME="${NAME%.gz}"
case "$NAME" in *.*) NAME="${NAME%.*}" ;; esac
BASE="$(dirname "$REF")/$NAME"

echo "Preparing reference: $REF"

echo "1. Creating FASTA index..."
samtools faidx "$REF"

echo "2. Creating sequence dictionary..."
samtools dict "$REF" -o "${BASE}.dict"

echo "3. Generating chromosome sizes..."
cut -f1,2 "${REF}.fai" > "${BASE}.chrom.sizes"

echo ""
echo "Files created:"
ls -la "${REF}.fai" "${BASE}.dict" "${BASE}.chrom.sizes"

echo ""
echo "Reference summary:"
echo "  Chromosomes: $(wc -l < "${REF}.fai")"
echo "  Total length: $(awk '{sum += $2} END {print sum}' "${REF}.fai") bp"
