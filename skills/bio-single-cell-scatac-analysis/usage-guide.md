# scATAC-seq Analysis - Usage Guide

## Overview

Single-cell ATAC-seq measures chromatin accessibility per cell, revealing cell-type-specific regulatory elements and TF activity. This skill centers on Signac (R/Seurat) with ArchR (R, on-disk, large data) and SnapATAC2 (Python, >1M cells) as alternatives. It treats the central epistemics of the data: a zero is ambiguous, the matrix is near-binary by sampling not biology, binarization is disfavored, gene activity is a weak proxy, and the peak set is circular with clustering.

## Prerequisites

See the Install section of `SKILL.md`.

## Quick Start

Tell your AI agent what you want to do:
- "Process my 10X scATAC fragments through QC, LSI, and clustering"
- "Diagnose which LSI components track sequencing depth and drop them"
- "Call consensus peaks per cell type and quantify a peak matrix"
- "Run chromVAR with a GC-matched background and rank TFs by z-score"
- "Detect homotypic and heterotypic doublets and combine them"
- "Pick a framework for 2 million cells"

## Example Prompts

### Data Loading and QC
> "Load my 10X filtered_peak_bc_matrix.h5 with the fragments file and add TSS enrichment and nucleosome signal"
> "Filter cells from the joint distribution of TSS enrichment and fragment count, not copied thresholds"

### Dimensionality Reduction
> "Run TF-IDF and SVD, then show DepthCor and drop only the components that track depth"
> "Use ArchR iterative LSI so rare populations are not lost in a single clustering pass"

### Peaks and Differential Accessibility
> "Call peaks per cluster on pseudobulk and merge to a fixed-width consensus set"
> "Find differentially accessible peaks with a logistic-regression test and fragment count as a latent variable"

### Motifs
> "Run chromVAR with GC-matched backgrounds and rank TFs by z-score"
> "This motif is enriched; confirm the actual TF with multiome RNA before claiming it drives the program"

### Doublets and Integration
> "Run AMULET for homotypic doublets and ArchR doublet scores for heterotypic, then combine"
> "Transfer cell-type labels from my scRNA-seq reference onto the ATAC cells"

## Related Skills

single-cell/multimodal-integration - joining the ATAC modality with RNA (Multiome WNN/MultiVI)
single-cell/preprocessing - shared QC and filtering concepts from scRNA-seq
single-cell/clustering - clustering and UMAP shared with scRNA-seq
single-cell/doublet-detection - doublet concepts and rate expectations
atac-seq/atac-peak-calling - bulk ATAC peak-calling background (MACS shift/extend)
atac-seq/motif-deviation - chromVAR deviation scoring in depth
chip-seq/motif-analysis - motif databases (JASPAR/cisBP) and enrichment testing
