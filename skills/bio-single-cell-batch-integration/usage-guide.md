# Batch Integration - Usage Guide

## Overview

Batch integration learns a shared representation for technical batches while preserving biological cell states. For the methods, safety boundaries, runnable commands, and evaluation criteria, use the matching sections of `SKILL.md`.

## Quick Start

Tell the agent what to do:
- "Integrate my samples and remove batch effects"
- "Pick an integration method for a large atlas with many batches"
- "Is my batch confounded with condition - should I integrate at all?"
- "Score my integration and check I did not over-correct"

## Example Prompts

### Integration
> "Merge my samples and run Harmony, then cluster on the corrected embedding"

> "Use scVI to integrate a large multi-batch atlas from raw counts"

> "Run Seurat v5 RPCA integration because the datasets only partly overlap"

### Method Choice
> "I have a few same-protocol samples - which method and should I even integrate?"

> "Some cells are labeled - would scANVI preserve biology better here?"

### Diagnosis
> "A rare population disappeared after integration - is this over-correction?"

> "My batch tracks my treatment - can integration separate them?"

### Scoring
> "Compute kBET, iLISI, and cell-type silhouette and tell me if biology was kept"

> "Benchmark Harmony vs scVI vs scANVI with scib-metrics and pick the most robust"

## Related Skills

- preprocessing - QC and normalization that must precede integration
- clustering - Cluster on the integrated embedding, not on raw PCA
- cell-annotation - Reference mapping and label transfer after integration
- single-cell/multimodal-integration - Joint analysis across modalities (distinct from batch integration)
- single-cell/differential-abundance - Test whether composition shifts across conditions after integration
- differential-expression/deseq2-basics - Pseudobulk DE on uncorrected counts per cell type
- data-visualization/dimensionality-reduction-plots - Before/after UMAP comparison figures
