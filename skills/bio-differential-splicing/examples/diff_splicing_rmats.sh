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
MIN_COVERAGE=10            # min inclusion + min skipping reads over all replicates (both groups)

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

# Filter significant events: FDR < 0.05, |deltaPSI| > 0.1, and >= MIN_COVERAGE reads (per-replicate minima summed)
echo ""
echo "Significant SE events (|deltaPSI| > 0.1, FDR < 0.05, coverage >= $MIN_COVERAGE):"
awk -F'\t' -v mincov="$MIN_COVERAGE" '
function minof(s,   a, n, i, m) { n = split(s, a, ","); m = a[1] + 0; for (i = 2; i <= n; i++) if (a[i] + 0 < m) m = a[i] + 0; return m }
NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; print; next }
$col["FDR"] != "NA" && $col["IncLevelDifference"] != "NA" {
    mi = minof($col["IJC_SAMPLE_1"]); t = minof($col["IJC_SAMPLE_2"]); if (t < mi) mi = t
    ms = minof($col["SJC_SAMPLE_1"]); t = minof($col["SJC_SAMPLE_2"]); if (t < ms) ms = t
    d = $col["IncLevelDifference"] + 0
    if ($col["FDR"] + 0 < 0.05 && (d > 0.1 || d < -0.1) && mi + ms >= mincov) print
}' "$OUTPUT_DIR/SE.MATS.JC.txt" | head -20
