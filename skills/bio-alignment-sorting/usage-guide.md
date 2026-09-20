# Alignment Sorting - Usage Guide

## Overview
Sort BAM files by coordinate or read name for downstream analysis tools that require specific sort orders. Commands, the sort order each downstream tool needs, and error fixes are in `SKILL.md`.

## Prerequisites
```bash
# samtools
conda install -c bioconda samtools

# pysam
pip install pysam
```

## Quick Start
Tell your AI agent what you want to do:
- "Sort my BAM file by coordinate"
- "Sort reads by name for duplicate marking"
- "Check the current sort order of my BAM file"
- "Sort directly from aligner output in a pipeline"

## Example Prompts

### Coordinate Sorting
> "Sort my BAM file by genomic position"

> "Sort the output from BWA alignment"

> "Verify my BAM is really coordinate-sorted, not just labelled that way"

### Name Sorting
> "Sort my BAM file by read name"

> "Sort by name so Picard MarkDuplicates accepts it"

> "Prepare BAM for the fixmate step"

> "Sort for extracting paired FASTQ files"

### Pipeline Operations
> "Align with BWA and sort in one command"

> "Run the complete duplicate marking workflow"

> "Re-sort an incorrectly sorted BAM file"

> "Merge these per-lane BAMs"

### Performance Optimization
> "Sort with 8 threads and 4GB memory per thread"

> "Use a fast SSD for temporary files during sorting"

## What the Agent Will Do

1. Check the current sort order of the input BAM (records, not just the header)
2. Select appropriate sort method (coordinate, `-n`/`-N` name, tag, template-coordinate)
3. Configure memory and thread settings for optimal performance
4. Execute the sort operation
5. Verify the output file was created successfully (`samtools quickcheck`)
6. Index the sorted file if coordinate-sorted

## Tips
- Keep the original unsorted file until you have verified the sorted output
- Related Skills: sam-bam-basics, alignment-indexing, duplicate-handling, alignment-filtering
