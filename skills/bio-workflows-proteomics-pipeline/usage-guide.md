# Proteomics Pipeline Usage Guide

## Overview

End-to-end workflow for label-free proteomics analysis, from MaxQuant/DIA-NN output to differential protein abundance, with TMT and SILAC variants. The agent follows `SKILL.md`; this guide is only for choosing the Skill.

## Quick Start

Tell your AI agent what you want to do:
- "Run the proteomics pipeline on my MaxQuant output"
- "Find differentially expressed proteins between conditions"
- "Process my DIA-NN results and run differential analysis"

## Example Prompts

### Basic Analysis
> "I have proteinGroups.txt from MaxQuant, run the full pipeline"

> "Normalize my proteomics data and find differential proteins"

### QC and Preprocessing
> "Check sample quality with PCA and correlation heatmap"

> "Handle missing values correctly by modeling the dropout rather than imputing"

### Differential Analysis
> "Run limma to find proteins changed between treatment and control"

> "Use MSstats for differential analysis with my peptide-level data"

> "Analyze my two-plex TMT experiment with the reference channel"

## Inputs and outputs

Inputs, the sample-annotation columns and the package install are in `SKILL.md` ("Inputs and Install").
`examples/proteomics_workflow.R` runs standalone on a simulated table and writes:

| File | Description |
|------|-------------|
| proteomics_results.csv | All proteins x contrasts with statistics |
| proteomics_results_raw_boxplot.pdf | Raw per-sample log2 distributions, inspected before normalization |
| proteomics_results_pca.pdf | Sample clustering |
| proteomics_results_volcano.pdf | Log2FC vs -log10(p-value) |
| proteomics_results_heatmap.pdf | Significant proteins |

## Typical Results

- 2000-5000 quantified proteins (cell lysate)
- 50-500 differential proteins (10%)
- Fold changes typically 1.5-4x

## Related Skills

See `SKILL.md` ("Related Skills").
