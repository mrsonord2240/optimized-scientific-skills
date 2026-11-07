#!/bin/bash
# Reference: CRISPResso2 2.2+, pandas 2.2+ | Verify API if version differs
# CRISPResso2 editing analysis example

set -e

# Configuration
AMPLICON="AATGTCCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAGGTCATGATCCCCTTCTGGAGCTCCCAACGGGCCGTGGTCTGGTTCATCATCTGTAAGAATGGCTTCAAGAGGCTCGGCTGTGGTT"
GUIDE="CTGCCCACAGGTGAGGAGGT"
OUTPUT_DIR="crispresso_results"

mkdir -p "$OUTPUT_DIR"

# Analyze control sample
echo "Analyzing control sample..."
CRISPResso \
    --fastq_r1 control_R1.fastq.gz \
    --fastq_r2 control_R2.fastq.gz \
    --amplicon_seq "$AMPLICON" \
    --guide_seq "$GUIDE" \
    --output_folder "$OUTPUT_DIR" \
    --name control \
    --quantification_window_size 10

# Analyze edited sample
echo "Analyzing edited sample..."
CRISPResso \
    --fastq_r1 edited_R1.fastq.gz \
    --fastq_r2 edited_R2.fastq.gz \
    --amplicon_seq "$AMPLICON" \
    --guide_seq "$GUIDE" \
    --output_folder "$OUTPUT_DIR" \
    --name edited \
    --quantification_window_size 10

# Compare samples
# CRISPRessoCompare takes the two output folders as positional arguments, not
# --crispresso_output_folder_1/_2 flags -- those flags do not exist and CRISPRessoCompare
# exits with "unrecognized arguments" if used (verified against CRISPResso2 2.3.4 --help).
echo "Comparing samples..."
CRISPRessoCompare \
    "$OUTPUT_DIR/CRISPResso_on_control" \
    "$OUTPUT_DIR/CRISPResso_on_edited" \
    --output_folder "$OUTPUT_DIR/comparison"

echo "Analysis complete!"
echo "Results in: $OUTPUT_DIR"
