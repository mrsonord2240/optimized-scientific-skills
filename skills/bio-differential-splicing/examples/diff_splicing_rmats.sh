#!/bin/bash
# Reference: rMATS-turbo 4.4.0, samtools 1.24 (checked 2026-09-20)
# Differential splicing analysis with rMATS-turbo
# Compares splicing between two conditions

set -e

# Configuration
GTF="annotation.gtf"
THREADS=8
OUTPUT_DIR="rmats_output"
TMP_DIR="rmats_tmp"
LIB_TYPE="fr-unstranded"   # fr-firststrand / fr-secondstrand for stranded libraries (check with RSeQC infer_experiment.py)
MIN_COVERAGE=10            # min inclusion + min skipping reads in at least half of each group's replicates

# Sample list files
# Each file has ONE line: comma-separated BAM paths for one condition
# Example: /path/to/sample1.bam,/path/to/sample2.bam,/path/to/sample3.bam
CONDITION1_BAMS="condition1_bams.txt"
CONDITION2_BAMS="condition2_bams.txt"

# Read type (-t) and read length come from the data, not from a hard-coded default:
# -t paired on single-end reads, or a --readLength that does not match the reads, exits 0 with an empty result table.
FIRST_BAM=$(head -n 1 "$CONDITION1_BAMS" | cut -d, -f1)
read -r READ_TYPE READ_LENGTH < <(samtools view "$FIRST_BAM" | head -n 10000 | \
    awk '{ if ($2 % 2) p++; if (length($10) > m) m = length($10) } END { print (p > 0 ? "paired" : "single"), m }')
echo "Detected from $FIRST_BAM: -t $READ_TYPE, --readLength $READ_LENGTH"

# Run rMATS-turbo
# Performs both quantification and differential testing
rmats.py \
    --b1 "$CONDITION1_BAMS" \
    --b2 "$CONDITION2_BAMS" \
    --gtf "$GTF" \
    -t "$READ_TYPE" \
    --readLength "$READ_LENGTH" \
    --variable-read-length \
    --libType "$LIB_TYPE" \
    --nthread "$THREADS" \
    --od "$OUTPUT_DIR" \
    --tmp "$TMP_DIR" \
    --cstat 0.01  # splicing-difference (dPSI) cutoff c for the null test |PSI1-PSI2|<=c (default 0.0001)

# rMATS exits 0 even when it produced nothing usable: check the table before trusting it
N_SE=$(( $(wc -l < "$OUTPUT_DIR/SE.MATS.JC.txt") - 1 ))
if [ "$N_SE" -lt 1 ] || ! head -n 1 "$OUTPUT_DIR/SE.MATS.JC.txt" | grep -q 'FDR'; then
    echo "ERROR: $OUTPUT_DIR/SE.MATS.JC.txt has no SE events or no FDR column (check -t, --readLength, --gtf contig names)" >&2
    exit 1
fi

# Output files:
# SE.MATS.JC.txt - Skipped exons (junction counts only)
# SE.MATS.JCEC.txt - Skipped exons (junction + exon body counts)
# A5SS, A3SS, MXE, RI - Alternative splice sites, mutually exclusive, retained introns

echo "rMATS analysis complete: $N_SE SE events. Results in $OUTPUT_DIR"
echo ""
echo "Key output columns:"
echo "  - IncLevelDifference: deltaPSI (positive = higher in condition1)"
echo "  - FDR: Benjamini-Hochberg corrected p-value"
echo "  - IJC/SJC: Inclusion/Skipping junction counts"

# Filter significant events: FDR < 0.05, |deltaPSI| > 0.1, and >= MIN_COVERAGE total reads in at least half of each group.
# Do not sum independent inclusion/skipping minima: they can come from different replicates and lose true calls as n grows.
echo ""
echo "Significant SE events (|deltaPSI| > 0.1, FDR < 0.05, coverage >= $MIN_COVERAGE):"
awk -F'\t' -v mincov="$MIN_COVERAGE" '
function half_or_more_covered(inc, skip, minimum,   a, b, n, i, covered) {
    n = split(inc, a, ","); if (split(skip, b, ",") != n) return 0
    for (i = 1; i <= n; i++) if (a[i] + b[i] >= minimum) covered++
    return covered >= int((n + 1) / 2)
}
NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; print; next }
$col["FDR"] != "NA" && $col["IncLevelDifference"] != "NA" {
    d = $col["IncLevelDifference"] + 0
    if ($col["FDR"] + 0 < 0.05 && (d > 0.1 || d < -0.1) &&
        half_or_more_covered($col["IJC_SAMPLE_1"], $col["SJC_SAMPLE_1"], mincov) &&
        half_or_more_covered($col["IJC_SAMPLE_2"], $col["SJC_SAMPLE_2"], mincov)) print
}' "$OUTPUT_DIR/SE.MATS.JC.txt" | head -20
