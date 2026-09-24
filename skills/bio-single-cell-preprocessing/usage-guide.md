# Single-Cell Preprocessing - Usage Guide

## Overview

This skill covers quality control, ambient-RNA handling, normalization, and feature selection for single-cell RNA-seq in Scanpy (Python) and Seurat (R). It is framed around the decisions that drive every downstream result: where to set QC thresholds without deleting real cell types, whether to remove ambient RNA and with which tool, which normalization to use, and whether to scale and regress out covariates.

## Example Prompts

### Quality Control
> "Compute QC metrics with mito, ribo, and hemoglobin fractions and show the distributions"

> "Filter cells using 5 MAD on counts/genes and a tissue-aware mitochondrial cap"

> "My tissue is cardiac muscle - set a mito threshold that does not delete cardiomyocytes"

### Ambient RNA
> "Run SoupX on my Cell Ranger output and report the contamination fraction"

> "Decide whether I need CellBender or DecontX for this snRNA-seq dataset"

### Normalization and Features
> "Normalize with shifted-log and explain whether I should use scran instead"

> "Select highly variable genes from raw counts with seurat_v3"

> "Should I scale and regress out total_counts before PCA?"

## Related Skills

- single-cell/data-io - load the raw matrix before preprocessing
- single-cell/doublet-detection - per-sample doublet calling around the QC step
- single-cell/clustering - PCA, neighbors, and clustering after preprocessing
- single-cell/batch-integration - correct batch effects instead of regressing them out
- single-cell/markers-annotation - find markers after clustering
- differential-expression/deseq2-basics - pseudobulk DE across samples (avoids single-cell pseudo-replication)
