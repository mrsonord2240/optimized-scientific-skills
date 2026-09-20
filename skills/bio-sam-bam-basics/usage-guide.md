# SAM/BAM/CRAM Basics - Usage Guide

## Overview
View, convert, and understand alignment files in SAM, BAM, and CRAM formats using samtools and pysam. Commands, tables (FLAG, CIGAR, MAPQ, tags), CRAM reference setup and failure modes live in `SKILL.md`; installation is in its "Version Compatibility" section.

## Quick Start
Tell your AI agent what you want to do:
- "View the first 10 reads in my BAM file"
- "Convert my SAM file to BAM format"
- "Count the total number of reads in sample.bam"
- "Extract reads from chromosome 1"

## Example Prompts

### Viewing Files
> "Show me the header of my BAM file"

> "View the first 20 alignments with the header included"

> "Count how many reads are in my BAM file"

### Format Conversion
> "Convert my SAM file to compressed BAM"

> "Convert my BAM to CRAM format using the reference genome"

> "Convert CRAM back to BAM for compatibility with older tools"

### Region Extraction
> "Extract all reads from chr1:1000000-2000000"

> "Get reads from multiple regions: chr1:1000-2000 and chr2:3000-4000"

> "Count reads in a specific genomic region"

### Understanding Alignments
> "Explain what FLAG 99 means in my SAM file"

> "Parse the CIGAR string 10M2I30M5D20M"

> "Show me the mapping quality distribution of my reads"

## What the Agent Will Do

1. Open the alignment file with appropriate mode (SAM/BAM/CRAM)
2. Parse the header to understand reference sequences and metadata
3. Iterate through alignments or fetch specific regions
4. Extract read properties (name, flag, position, quality, sequence)
5. Apply any requested filters or transformations
6. Output results in the requested format

## Related Skills
See the "Related Skills" list at the end of `SKILL.md` (indexing, sorting, filtering, validation, statistics, reference operations).
