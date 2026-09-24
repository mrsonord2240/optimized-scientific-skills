# Proteomics QC - Usage Guide

## Overview
Bottom-up proteomics QC as a three-level funnel (instrument, identification/run, experiment/quantitative). The agent inspects raw signal and removes contaminants BEFORE normalizing, because median normalization erases loading failures. Metrics, thresholds, code and failure modes are in `SKILL.md`.

## Prerequisites
```bash
pip install numpy pandas scipy scikit-learn matplotlib seaborn
# R packages: install.packages('PTXQC')
# BiocManager::install(c('limma', 'MSstatsTMT'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Plot raw per-sample boxplots and ID counts before I normalize, and flag loading failures"
- "Strip MaxQuant contaminant and decoy rows before log-transform and normalization"
- "Compute replicate correlation on log2 and CV on the linear scale, then run PCA colored by batch"
- "Diagnose whether my missing values are MNAR or MCAR before I pick an imputer"

## Example Prompts

### Inspect Before Normalizing
> "Show raw, un-normalized boxplots per sample with ID counts and total signal, and flag any sample shifted more than 2-3x below its group"

> "Remove rows flagged Potential contaminant, Reverse, or Only identified by site from my MaxQuant proteinGroups before I normalize"

> "What fraction of summed intensity is keratin and trypsin, and is it higher in my low-input samples?"

### Reproducibility
> "Calculate within-group Pearson correlation on log2 intensities and check whether any sample matches a different group better (possible swap)"

> "Compute median CV per condition on the linear scale, not on log data"

> "Are my technical replicates above r 0.98 and is the biological CV in the 20-40% range?"

### Missing Values
> "Plot present-fraction versus abundance to decide if missingness is MNAR or MCAR"

> "Filter to proteins valid in at least 70% of replicates in one condition before imputing"

> "Which imputation method matches my missingness mechanism, and why would the wrong one corrupt my results?"

### Batch and Outliers
> "Run PCA on the normalized survivors and test whether PC1 or PC2 associates with batch rather than condition"

> "Is batch the dominant axis of variance, and should I correct it before differential testing?"

### TMT and DIA
> "Check retention-time fit and peak width per run in my DIA-NN report"

> "Check TMT channel-loading balance on raw reporter intensities and flag any channel deviating more than 2x"

> "How many protein groups pass at 1% global q-value, and why is precursor q not enough?"

## Related Skills
- data-import - Load search-engine output and intensity matrices before QC
- quantification - Normalization and imputation mechanics that QC mandates running AFTER inspection
- differential-abundance - The moderated statistical test QC gates
- dia-analysis - DIA q-value/FDR internals behind the protein-count QC
- data-visualization/dimensionality-reduction-plots - PCA/MDS projection plotting
- workflows/proteomics-pipeline - End-to-end pipeline placing QC before differential testing
