# Dimensionality-Reduction Plots - Usage Guide

## Overview

PCA, t-SNE, UMAP, and PHATE are the four standard methods for projecting high-dimensional omics data to 2D. Each preserves a different property - variance (PCA), local neighborhoods (t-SNE), local + partial global (UMAP), continuous transitions (PHATE). The Chari-Pachter 2023 critique established that 2D embeddings lose >95% of high-dimensional geometry, so embeddings communicate "these points are similar locally" and nothing more. Hyperparameters, random seeds, and explicit interpretation limits matter for reproducibility.

## Prerequisites

```bash
pip install scanpy umap-learn openTSNE phate scikit-learn matplotlib scikit-misc igraph
```

The shipped `examples/embedding_phd.py` expects an `.h5ad` with raw, non-negative integer counts in `X`, at least 100 cells and 2,000 genes, and `obs['condition']`. `obs['pseudotime']` is optional. It uses Scanpy's igraph Leiden flavor, so it does not require the separate `leidenalg` Python package.

```r
install.packages(c('Rtsne', 'uwot', 'PCAtools', 'phateR'))
BiocManager::install(c('PCAtools', 'DESeq2'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Make a PCA plot of bulk RNA-seq for sample QC, colored by condition and shaped by batch"
- "Compute UMAP from scanpy AnnData using n_neighbors=30, min_dist=0.3, random_state=42"
- "Run openTSNE with an explicit PCA init, perplexity 30, auto learning rate, and random_state 42"
- "Plot PHATE for trajectory display instead of UMAP"
- "Annotate PCA axes with variance explained percentages"
- "Show loadings as arrows on PC1 vs PC2"

## Example Prompts

### Bulk PCA for sample QC

> "PCA on vst-normalized bulk RNA-seq counts. Plot PC1 vs PC2 colored by condition, shaped by batch. Axis labels must include variance explained. Add screeplot of PC1-10."

### Single-cell UMAP

> "From scanpy AnnData, compute PCA(50), neighbors(30, n_pcs=50), UMAP(min_dist=0.3, random_state=42). Plot colored by Leiden cluster with on-data labels."

### Kobak-Berens t-SNE

> "Run openTSNE with PCA initialization, perplexity=30, learning_rate='auto', and random_state=42. Plot the embedding colored by cluster assignment."

### PHATE for trajectory

> "Use PHATE instead of UMAP for displaying developmental data - preserves continuous transitions."

### Method comparison

> "Run PCA, UMAP, t-SNE, and PHATE on the same matrix. Display side-by-side. Annotate each panel with its preservation property."

## What the Agent Will Do

1. Decide method: PCA for variance-explained QC; UMAP for cluster overview; t-SNE if cluster boundaries are critical; PHATE for continuous trajectories.
2. Pre-process input: normalize library size, log-transform, then scale; compute PCA(50) before t-SNE/UMAP for single-cell. Run `seurat_v3` HVG selection on raw counts, before normalization.
3. Set explicit hyperparameters: perplexity / n_neighbors / min_dist / random_state.
4. Fit the projection with a fixed seed for reproducibility.
5. Plot with axis labels: PCA shows variance %; UMAP/t-SNE/PHATE label only "UMAP1 / UMAP2" (no units).
6. Color by categorical (CVD-safe palette) or continuous (perceptually-uniform colormap).
7. Annotate with cluster labels on-plot OR via legend depending on cluster count.
8. State in caption: hyperparameters used; embedding's interpretation limit.

## Implementation Notes

All thresholds, method-specific recipes, failure modes, seed semantics, output-path guidance, and interpretation limits live in `SKILL.md`. Its **Reference Files** section links the runnable comparison and the neighborhood/batch validation procedure; use those instead of duplicating the rules here.

## Related Skills

- single-cell/preprocessing - PCA / neighbor graph before embedding
- single-cell/clustering - Leiden / Louvain assignments to color UMAP
- single-cell/trajectory-inference - Pseudotime / RNA velocity for trajectory claims
- data-visualization/color-palettes - Categorical and perceptual palettes
- data-visualization/distribution-plots - Per-cluster gene-expression follow-up
