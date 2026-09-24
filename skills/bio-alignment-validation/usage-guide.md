# Alignment Validation - Usage Guide

## Overview
Validate alignment quality with insert size distribution, proper pairing rates, GC bias, strand balance, and other post-alignment metrics before downstream analysis.

## Prerequisites
samtools, Picard and pysam; install line in SKILL.md "Version Compatibility".

## Quick Start
Tell your AI agent what you want to do:
- "Validate the quality of my aligned BAM file"
- "Check insert size distribution for sample.bam"
- "Generate a comprehensive QC report for my alignment"
- "Check for GC bias in my sequencing data"

## Example Prompts

### Basic Validation
> "Run alignment validation on sample.bam and tell me if it passes QC"

> "Check the mapping rate and proper pairing percentage for my BAM file"

> "Generate a flagstat report and interpret the results"

### Insert Size Analysis
> "Calculate insert size distribution for sample.bam and plot a histogram"

> "Check if the insert size matches my library prep protocol (expected 350bp)"

> "Compare insert sizes across multiple samples to check for consistency"

### GC Bias Detection
> "Check for GC bias in my WGS data using Picard"

> "Generate a GC bias plot and tell me if correction is needed"

> "Compare GC coverage across my samples to identify problematic ones"

### Strand Balance
> "Calculate the forward fraction F/(F+R) for sample.bam"

> "Check strand balance per chromosome to identify bias"

> "Validate that the forward fraction F/(F+R) is near 0.5"

### Comprehensive QC
> "Run a complete alignment validation pipeline with all metrics"

> "Generate QC report with mapping rate, pairing, insert size, and strand balance"

> "Identify any quality issues in my alignment before variant calling"

## What the Agent Will Do

1. Check the BAM is intact (`samtools quickcheck`; the example validators need no index)
2. Calculate mapping statistics using samtools flagstat
3. Extract insert size distribution from properly paired reads
4. Compute strand balance (forward fraction F/(F+R))
5. Check mapping quality distribution
6. Generate per-chromosome coverage statistics
7. Optionally run Picard metrics for GC bias and alignment summary
8. Compare metrics against quality thresholds
9. Generate a summary report with pass/fail status
10. Recommend corrective actions for any failing metrics

## Where the Details Live
Thresholds, commands, interpretation of each metric and what to check when one fails are in SKILL.md (Quality Thresholds Summary, When a Metric Fails); related Skills are listed at its end. Ready-to-run validators: `examples/validate_alignment.sh` and `examples/validate_alignment.py` (exit 0 pass/warn, 1 fail, 2 unreadable or empty).
