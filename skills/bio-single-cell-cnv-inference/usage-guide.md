# Copy-Number Inference from Single-Cell RNA-seq - Usage Guide

## Overview

This skill infers large-scale copy-number alterations from tumor single-cell or single-nucleus RNA-seq by treating averaged expression over genomic windows as a proxy for DNA copy number. It separates malignant from normal cells and calls subclones at chromosome-arm / large-segment (~5 Mb) resolution. It covers reference-based expression smoothing (inferCNV), reference-free segmentation (copyKAT, SCEVAN), and haplotype-aware allele-plus-expression inference (Numbat). It is distinct from DNA/WES-based copy-number (copy-number/cnvkit-analysis), which measures DNA read depth and resolves focal events.

## Quick Start

Tell your AI agent what you want to do:
- "Which cells in my tumor scRNA-seq are malignant?"
- "Infer chromosome-arm CNVs from my tumor expression matrix"
- "Run inferCNV using my T cells and myeloid cells as the normal reference"
- "Call tumor subclones with allele-aware Numbat"
- "Pick a CNV-inference method when I have no annotated normal cells"

## Example Prompts

### Malignant vs normal
> "Separate malignant from normal cells in this tumor sample using a CNV-inference method"
> "I annotated immune and stromal cells; use them as the inferCNV reference and call which cells are aneuploid"

### Method choice
> "I have no normal cells annotated - which reference-free CNV caller should I use?"
> "Should I use inferCNV, copyKAT, or Numbat for this dataset, and why?"

### Subclones and alleles
> "Run Numbat to resolve subclones and copy-neutral LOH from my BAM and counts"
> "My copyKAT subclones look unstable - how do I confirm them?"

### Interpretation and pitfalls
> "My immune cells are being called aneuploid - what went wrong with my reference?"
> "Run CNV inference per patient instead of on my integrated cross-patient object"

## Details

Install commands, thresholds, failure modes and interpretation live in `SKILL.md` (Install, Threshold and parameter reference, Common Errors).

## Related Skills

- single-cell/preprocessing - QC and normalization that precede CNV inference
- single-cell/clustering - Provides the cell groups and the malignant-vs-normal clustering the CNV call refines
- single-cell/cell-annotation - Identifies the non-malignant lineages used as the normal reference
- single-cell/batch-integration - Why CNV inference must run per patient before any cross-patient integration
- copy-number/cnvkit-analysis - DNA/WES-based copy-number that resolves focal events, the orthogonal contrast to this expression proxy
