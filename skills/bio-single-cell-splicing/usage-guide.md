# Single-Cell Splicing - Usage Guide

## Overview
Analyze alternative splicing at single-cell resolution. The first decision is library chemistry, not tool: 10X 3' cannot support transcriptome-wide splicing, plate-based full-length and single-cell long-read data can. Tools covered: MARVEL, BRIE2, scQuint, SpliZ, Psix, Sierra (APA, not splicing) and pseudobulk leafcutter. The chemistry gate, thresholds, and failure-mode index are in `SKILL.md`; runnable tool workflows are in `references/`.

## Example Prompts

### Chemistry Audit
> "I have 10X Chromium 3' v3 data; can I do splicing analysis? What can I do instead?"

### Plate-Based Workflow
> "Run MARVEL on Smart-seq2 BAMs to quantify SE/A5SS/A3SS PSI, classify modality, and test differential splicing between cell types."

> "Use BRIE2 with sequence-feature prior to estimate per-cell PSI for cassette exons in low-coverage Smart-seq3 data."

### Discovery
> "Run SpliZ to find genes with cell-state-associated splicing without using an event database."

### Trajectory
> "Use Psix to detect regulated AS along a developmental pseudotime trajectory, robust to dropout."

### APA (10X 3')
> "Use Sierra to peak-call 3' ends and detect alternative polyadenylation across cell types."

### Pseudobulk
> "Aggregate cells by donor and cluster and run leafcutter on the pseudobulk junction counts for differential splicing between cell types."

### Long-Read Single-Cell
> "I have MAS-Iso-seq + 10X 5' data; demultiplex barcodes with skera/lima, then run FLAMES for full-length isoform quantification per cell."

## Related Skills

- single-cell/preprocessing - QC and normalization
- single-cell/clustering - Cell type annotation prerequisite
- single-cell/data-io - h5ad / Seurat I/O
- splicing-quantification - Bulk RNA-seq comparison
- long-read-splicing - Full-isoform analysis from MAS-Iso-seq
