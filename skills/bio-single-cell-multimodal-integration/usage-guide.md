# Multimodal Integration - Usage Guide

## Overview

Multimodal single-cell assays measure several biological layers per cell (RNA + surface protein in CITE-seq, RNA + chromatin accessibility in 10x Multiome) or pair independent single-modality datasets (unpaired scRNA + scATAC, or a mosaic mix of both). This skill classifies the integration task by anchor structure, denoises each modality natively, and selects a joint method (WNN, totalVI, MultiVI, MOFA+, GLUE, Seurat v5 bridge) with explicit failure modes. Install commands are in SKILL.md's Prerequisites section.

## Quick Start

Tell your AI agent what you want to do:
- "Classify my multimodal task as paired or unpaired and pick a method"
- "Denoise my CITE-seq ADT with DSB before clustering"
- "Run WNN on my CITE-seq data and show me the modality weights"
- "Integrate independent scRNA and scATAC that have no shared cells"
- "Decide between totalVI and WNN for my CITE-seq experiment"

## Example Prompts

### CITE-seq (RNA + Protein)
> "My ADT clusters look like background smears; denoise with DSB using empty droplets, then run WNN"
> "Train totalVI on my CITE-seq MuData and give me denoised protein plus foreground probabilities"
> "WNN clustering seems driven by CD markers only; show the per-cell modality weight distribution and test stability if ADT is down-weighted"

### Multiome (RNA + ATAC, same cell)
> "Process RNA with PCA and ATAC with TF-IDF/LSI, drop depth-correlated LSI components, then join with WNN"
> "I merged two multiome runs and see batch structure; re-quantify against a unified peak set"

### Unpaired / Diagonal / Mosaic
> "Align my independent scRNA and scATAC with GLUE using a peak-near-gene guidance graph"
> "Map my scATAC query onto an scRNA reference using a multiome bridge dataset"
> "Integrate a mosaic design where one batch has RNA+ATAC and another has RNA only"

## What the Agent Will Do

1. Classify the task by anchor structure (vertical/paired, diagonal/unpaired, mosaic) to narrow the method class.
2. Run per-modality QC, denoise ADT, reduce each modality natively, then run the joint method that matches the anchor structure.
3. Report per-cell modality weights and flag any imputed modality as inference, not measurement.

See SKILL.md's Governing Principle, Classify the Task table, Method Decision Tables, and Common Errors table for the full anchor-structure logic, method tradeoffs, ADT normalization choice (DSB vs CLR), and documented failure modes -- that is the single source for all of it.

## Related Skills

single-cell/scatac-analysis - ATAC QC, TF-IDF/LSI, gene-activity caveats for the Multiome ATAC half
single-cell/preprocessing - per-modality RNA QC and normalization before integration
single-cell/clustering - clustering and UMAP on the joint graph
single-cell/batch-integration - horizontal (same-modality, cross-sample) correction
single-cell/markers-annotation - marker-based interpretation of joint clusters
atac-seq/motif-deviation - chromVAR TF activity on the Multiome ATAC modality
pathway-analysis/go-enrichment - functional interpretation of modality-specific factors
