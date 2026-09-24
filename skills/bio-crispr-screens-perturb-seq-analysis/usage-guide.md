# Perturb-Seq Analysis - Usage Guide

## Overview

Decision-grade analysis of single-cell pooled CRISPR screens. Covers Perturb-seq, CROP-seq, Perturb-CITE-seq, ECCITE-seq, Perturb-ATAC, and Multiome variants; MOI considerations; sgRNA assignment; escaper filtering with Mixscape (Papalexi 2021); calibrated low-MOI DE with SCEPTRE (Barry 2024); Pertpy unified framework; and genome-wide screens (Replogle 2022).

## Prerequisites

```bash
# Pertpy (Python; integrates Mixscape + SCEPTRE wrappers + DE methods)
pip install pertpy scanpy anndata muon

# SCEPTRE (R; for calibrated low-MOI DE)
R -e "devtools::install_github('katsevich-lab/sceptre')"   # sceptre is not on CRAN

# Seurat (R; for Mixscape native + multiome)
R -e "install.packages('Seurat')"
```

Required inputs:
- AnnData / Seurat object with scRNA-seq counts (cells x genes)
- Per-cell sgRNA assignment (from direct-capture or CROP-seq library)
- Sample / channel metadata (batch covariate)
- Non-targeting control sgRNA identifier

## Quick Start

Tell the AI agent what to do:
- "Analyze my CROP-seq experiment: sgRNA assignment, Mixscape escaper filter, SCEPTRE DE per perturbation, downstream pathway analysis"
- "Choose between direct-capture Perturb-seq and CROP-seq architecture for a planned 1000-pert screen in iPSC-derived neurons"
- "Scale up: design genome-wide Perturb-seq following the Replogle 2022 protocol (2.5M cells, ~9,866 expressed genes in K562, CRISPRi)"
- "Diagnose: why does my Mixscape filter out 80% of perturbed cells as escapers?"
- "Run SCEPTRE on my low-MOI Perturb-seq data and compare FDR calibration vs MAST"

## Example Prompts

### Architecture Selection

> "I'm running a 1,500-perturbation screen in primary T cells. Choose between Dixit Perturb-seq (direct capture) vs CROP-seq vs Perturb-CITE-seq. Required: scRNA + sgRNA detection from same library. Recommendation depends on cost vs sgRNA assignment rate."

> "Genome-wide essentiality screen in K562. Replicate Replogle 2022 design: 2,057 essential-gene perturbations via CRISPRi, 10X 3' direct capture, a median >100 cells/pert as screened (budget 500-1,000 if per-perturbation DE power is the goal). Estimate cost and channel count."

### sgRNA Assignment + Filtering

> "Assign sgRNAs per cell using threshold of 10 reads. Compute multiplet rate (cells with 2+ sgRNAs above threshold); flag for doublet filter."

> "Apply Mixscape to filter escapers. For each perturbation, compute perturbation signature = cell_expression - mean(K=20 nearest NTC cells). Classify KO vs NP cells. Keep only KO cells for DE."

### DE Analysis

> "Run SCEPTRE on KO-filtered cells. NB GLM with technical covariates (n_genes, n_umi, channel). Permutation FDR (1000 iterations). Output per-gene-per-pert log-fold-change + FDR."

> "Compare SCEPTRE vs MAST on the same Perturb-seq data. Show FDR calibration via permutation: expect MAST inflation; SCEPTRE calibrated."

### Genome-Wide Scale

> "Design genome-wide Perturb-seq for 19,000 protein-coding genes using one dual-sgRNA CRISPRi element per gene. Compute cells needed at 500/pert. Distribute across 10X channels."

> "For my 2-million-cell genome-scale dataset, run per-pert SCEPTRE; output a 2057-pert x 19000-gene matrix of log-FCs. Cluster perturbations by their gene-effect profiles to identify functional modules."

### Diagnostics

> "Mixscape filtered 78% of perturbed cells as escapers. Diagnose: weak phenotype (perturbation insufficient), wrong K parameter, or true escaper rate is high (Cas9 expression heterogeneous)?"

> "sgRNA assignment rate is 65% of cells (35% unassigned). Diagnose: under-loaded sgRNA library, wrong read threshold, or library-prep architecture mismatch?"

> "MAST DE called 8,000 significant genes per perturbation. Verify against SCEPTRE; SCEPTRE will likely call 500-1,500 -- the difference is MAST's uncalibrated FDR."

### Multi-Omic

> "I have 10X Multiome data with sgRNA capture (RNA + ATAC). Use muon for joint analysis; identify perturbation-specific chromatin + RNA changes."

> "Perturb-CITE-seq: integrate sgRNA + scRNA + surface ADT. Identify perturbations that change cell-surface phenotype."

## What the Agent Will Do

1. Identify experimental architecture (direct-capture vs CROP-seq vs Perturb-CITE vs Multiome)
2. Verify sgRNA library prep matches architecture
3. Standard scRNA QC: gene counts, UMI counts, mitochondrial %, doublet detection (Scrublet)
4. sgRNA assignment per cell; compute assignment rate
5. Flag multiplets (cells with 2+ sgRNAs); decide to filter or analyze as combinatorial
6. Standard normalization (scanpy: total + log1p; or scran)
7. Apply Mixscape escaper filtering with NTC controls and K=20 nearest neighbors
8. Verify KO retention rate (strongly guide-dependent -- see SKILL.md's Quantitative Thresholds)
9. Per-perturbation DE via SCEPTRE (low-MOI variant if applicable) with covariates (n_genes, n_umi, channel)
10. Permutation FDR (see SKILL.md's Quantitative Thresholds)
11. Aggregate per-perturbation signatures; pathway enrichment
12. For genome-scale: cluster perturbations by effect profiles
13. Output: per-pert DE tables, perturbation cluster heatmap, pathway analysis

## Tips

- For combinatorial Perturb-seq (intentionally high MOI), pair guide-pairs cassettes and analyze as combinatorial (see crispr-screens/combinatorial-screens).
- For everything else -- architecture choice, MOI/assignment thresholds, Mixscape, SCEPTRE vs MAST, cells-per-perturbation targets -- see SKILL.md's Experimental Architecture Comparison, MOI and sgRNA Assignment, Escaper Cell Filtering, Failure Modes, and Quantitative Thresholds sections.

## Validation Checklist

See SKILL.md's Quantitative Thresholds and Common Errors tables.

## Related Skills

- crispr-screens/library-design - Direct-capture vs CROP-seq library design
- crispr-screens/screen-qc - sgRNA assignment as QC
- crispr-screens/mageck-analysis - Pseudobulk alternative
- crispr-screens/hit-calling - Pseudo-bulk hit calling
- single-cell/preprocessing - scRNA-seq preprocessing
- single-cell/clustering - Post-DE clustering of perturbations
- single-cell/multimodal-integration - Multiome Perturb-seq
- single-cell/perturb-seq - General single-cell screen tools
- pathway-analysis/go-enrichment - Pathway analysis
