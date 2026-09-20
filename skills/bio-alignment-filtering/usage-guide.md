# Alignment Filtering - Usage Guide

## Overview
Filter alignments by mapping quality, flags, regions, and other criteria to prepare clean datasets for downstream analysis.

## Prerequisites
```bash
# samtools
conda install -c bioconda samtools

# pysam
pip install pysam
```

## Quick Start
Tell your AI agent what you want to do:
- "Remove unmapped reads from my BAM file"
- "Filter for high-quality mapped reads (MAPQ >= 30)"
- "Extract reads from specific genomic regions"
- "Remove duplicates and secondary alignments"

## Example Prompts

### Quality Filtering
> "Keep only reads with mapping quality 30 or higher"

> "Remove all unmapped reads from my BAM file"

> "Filter for primary alignments only (no secondary or supplementary)"

### Flag-Based Filtering
> "Remove duplicate reads from my BAM file"

> "Keep only properly paired reads"

> "Extract first-in-pair reads only"

### Region Filtering
> "Extract reads from chr1:1000000-2000000"

> "Get all reads overlapping my target BED file"

> "Extract reads from multiple regions"

### Subsampling
> "Subsample my BAM to 10% of reads"

> "Downsample to approximately 1 million reads"

> "Create a reproducible subset for testing"

## What the Agent Will Do

1. Analyze the current BAM file to understand read composition (and, before removing duplicates, whether they are marked at all)
2. Apply requested filters using appropriate flags and quality thresholds
3. Count reads before and after filtering
4. Write filtered output to a new BAM file
5. Verify the output file and report filtering statistics

FLAG tables, MAPQ thresholds per aligner, assay-specific filters, region and subsampling recipes, and the pysam versions are in `SKILL.md`; `examples/filter_bam.py` is a ready-made pysam command-line filter.

## Tips
- Tell the agent which aligner made the BAM: "drop multi-mapped reads" needs a different `-q` for BWA, Bowtie2, HISAT2 and STAR
- Tell the agent which caller comes next: SV callers need supplementary alignments that the standard filter removes
