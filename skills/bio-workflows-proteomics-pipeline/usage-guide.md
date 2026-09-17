# Proteomics Pipeline Usage Guide

## Overview

End-to-end workflow for label-free proteomics analysis from MaxQuant/DIA-NN output to differential protein abundance.

## Prerequisites

```r
BiocManager::install(c('limma', 'DEqMS', 'proDA', 'MSstats', 'MSstatsTMT', 'ashr'))
install.packages(c('pheatmap', 'ggplot2', 'arrow', 'iq'))
```

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

## Pipeline Stages

### 1. Data Import
- Load proteinGroups.txt (MaxQuant) or report.parquet (DIA-NN 1.9+)
- Filter contaminants, reverse, and only-identified-by-site BEFORE normalizing
- Extract intensity columns (LFQ intensity, not raw Intensity, for between-sample work) -- but inspect the raw `Intensity` columns for failed loads: MaxLFQ has already renormalized the LFQ columns, so a low injection is largely invisible there

### 2. Transformation
- Replace 0 with NA
- Log2 transform
- Median centering normalization

### 3. Completeness Filtering
- Keep proteins valid in >= ~60% of replicates in at least one condition
- Filter before any missing-value handling, so only shallow gaps remain

### 4. Missing-Value Handling
- Preferred: MODEL the left-censored MNAR dropout (proDA / msqrob2 / MSstats-AFT); no imputation
- Fallback only: downshift imputation, which manufactures false positives for on/off proteins near the detection limit
- Report a protein fully missing in one group as "undetected", not as a fold change

### 5. Quality Control
- PCA: Check replicate clustering
- Correlation heatmap: Sample similarity
- Missing value patterns: Random or systematic

### 6. Differential Analysis
- limma: empirical-Bayes moderated t-test (trend=TRUE, robust=TRUE); treat() for a minimum fold change
- DEqMS: count-aware moderation when quantification depth varies
- proDA: probabilistic dropout model (no imputation)
- MSstats / MSstatsTMT: feature-level mixed models, plus the IRS bridge for multi-batch TMT

### 7. Output
- Differential proteins table
- Volcano plot
- Heatmap of significant proteins

## Input Requirements

### MaxQuant Output
```
proteinGroups.txt  # Protein-level quantification
evidence.txt       # Peptide-level (for MSstats)
annotation.csv     # Sample metadata
```

### Sample Annotation
`batch` is required whenever samples were acquired in more than one run/day/plex -- the pipeline's
design branch keys on that column and puts batch in the model as a covariate; without it the batch
effect stays in the residual. `condition` may have more than two levels (dose series, time course);
every non-reference level is contrasted against the first level.

```csv
sample,condition,replicate,batch
Sample1,Control,1,B1
Sample2,Control,2,B2
Sample3,Treatment,1,B1
Sample4,Treatment,2,B2
```

## Expected Outputs

| File | Description |
|------|-------------|
| proteomics_results.csv | All proteins with statistics |
| proteomics_results_volcano.pdf | Log2FC vs -log10(p-value) |
| proteomics_results_raw_boxplot.pdf | Raw per-sample log2 distributions, inspected before normalization |
| proteomics_results_pca.pdf | Sample clustering |
| proteomics_results_heatmap.pdf | Significant proteins |

## Typical Results

- 2000-5000 quantified proteins (cell lysate)
- 50-500 differential proteins (10%)
- Fold changes typically 1.5-4x

## Tips

- **Missing values**: model the MNAR dropout (proDA/msqrob2/MSstats-AFT) rather than impute; downshift imputation manufactures systematic false positives
- **Normalization**: median centering assumes most proteins are unchanged AND that the changes are roughly symmetric up/down; on a one-sided design it shifts the whole unchanged proteome the other way and manufactures one-directional hits. Never normalize an AP-MS/enrichment pulldown this way; normalize on a spike-in or unchanged-protein set instead
- **More than two conditions**: supported -- contrasts are built from the condition levels against the reference (first) level, and adjusted across all contrasts with `decideTests(method = 'global')`
- **Completeness**: filter on per-group completeness before any missing-value handling, not a blanket >50% rule
- **Replicates**: minimum 3 biological replicates per condition
- **Contaminants**: filter MaxQuant contaminants, reverse, and only-identified-by-site BEFORE log + normalize
