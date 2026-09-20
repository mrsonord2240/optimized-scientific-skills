# Duplicate Handling - Usage Guide

## Overview
Mark and remove PCR and optical duplicates from alignment files to prevent bias in downstream variant calling and peak detection. All commands, thresholds, the platform `-d` table, error messages and the assay decision table live in `SKILL.md`.

## Prerequisites
```bash
# samtools
conda install -c bioconda samtools

# pysam
pip install pysam
```
UMI libraries also need `umi_tools` and `fgbio`; see "UMI-Aware Deduplication" in `SKILL.md`.

## Quick Start
Tell your AI agent what you want to do:
- "Mark duplicates in my BAM file"
- "Remove PCR duplicates from my alignment"
- "Calculate the duplicate rate for my sample"
- "Run the full duplicate marking workflow"

## Example Prompts

### Marking Duplicates
> "Mark duplicates in sample.bam and keep them flagged"

> "Remove duplicates completely from my BAM file"

> "Mark duplicates and generate statistics"

### Workflow Operations
> "Run the complete fixmate and markdup pipeline"

> "Prepare my BAM for duplicate marking with fixmate"

> "Mark duplicates with optical duplicate detection for NovaSeq data"

### UMI and Special Assays
> "My capture BAM has UMIs in the RX tag; deduplicate it and call consensus reads"

> "Should I mark duplicates in my RNA-seq / amplicon BAM?"

### Quality Assessment
> "Calculate the duplicate rate for my sample"

> "Check if my duplicate rate is acceptable for WGS"

> "Compare duplicate rates across my samples"

## What the Agent Will Do

0. Confirm the assay; if standard markdup is wrong for it (RNA-seq, amplicon, UMI, scRNA), stop and use the tool `SKILL.md` names
1. Name-sort (or collate), add mate information with `fixmate -m`, re-sort by coordinate
2. Mark (or remove) duplicates with markdup, with the platform's optical distance
3. Index the final BAM and check the record count
4. Report the duplicate rate and assess quality

## Related Skills
See "Related Skills" in `SKILL.md`.
