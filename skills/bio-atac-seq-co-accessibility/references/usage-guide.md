# Co-accessibility (cis-Regulatory Linkage) - Usage Guide

## Overview

Infer enhancer-promoter and enhancer-enhancer cis-regulatory connections from scATAC-seq using cell-to-cell variability in joint peak accessibility. Cicero, ArchR getCoAccessibility, and SCENIC+ produce peak-pair connection scores; thresholding generates enhancer-gene candidates. Co-accessibility is a hypothesis generator -- it overlaps Hi-C/Micro-C 3D contacts only partly (dataset dependent) and is not equivalent. Validate against orthogonal data (Hi-C, CRISPRi-FlowFISH).

## Prerequisites

```r
BiocManager::install(c('GenomicInteractions', 'Gviz'))
# monocle3 and the monocle3-compatible Cicero (1.3.x, DESCRIPTION Depends: monocle3) are GitHub-only;
# the Bioconductor `cicero` targets monocle2 and lacks new_cell_data_set/preprocess_cds
remotes::install_github('cole-trapnell-lab/monocle3')
remotes::install_github('cole-trapnell-lab/cicero-release', ref='monocle3')
remotes::install_github('GreenleafLab/ArchR', ref='master')
```

```bash
# SCENIC+ (1.0a2 alpha, Python <=3.11.8) for TF-driven networks (Multiome only); source-only install, pulls in pycisTopic. Not executed in this Skill's audit
git clone https://github.com/aertslab/scenicplus && cd scenicplus && pip install .
```

Inputs: scATAC peak-cell matrix (from Signac, ArchR, or SnapATAC2 preprocessing) and a UMAP/dim-reduction for metacell construction. SCENIC+ additionally requires paired RNA AnnData (Multiome).

## Quick Start

Tell your AI agent what you want to do:
- "Run Cicero on a Signac scATAC-seq object to get peak-pair co-accessibility"
- "Use ArchR getCoAccessibility for an ArchR project; connection cutoff 0.5"
- "Find enhancer-gene candidates by overlapping connections with promoters (TSS +/- 2 kb)"
- "Compare Cicero connections against published Hi-C loops to estimate concordance"
- "Run SCENIC+ on Multiome data to infer TF -> enhancer -> gene networks"

## Example Prompts

### Standard Cicero
> "Convert my Signac scATAC peak matrix to a CDS, reduce dimensions with LSI + UMAP, build a Cicero CDS via metacells (k=50), then run_cicero with `window=500000` for cis connections. Filter to coaccess > 0.25 as high-confidence."

### Per-Cluster Co-accessibility
> "Run Cicero separately on each Seurat cluster to capture cell-type-specific connections; merge results with cluster annotation."

### ArchR Wrapped Cicero
> "In my ArchR project, addCoAccessibility on PeakMatrix with k=100, knnIteration=500, maxDist=250000 (250 kb cis). Use corCutOff=0.5 for high-confidence."

### Enhancer-Gene Linking
> "Use the strong co-accessibility connections; map one anchor to gene promoters (TSS +/- 2 kb); the other anchor is the candidate enhancer. Output enhancer-gene pairs."

### Hi-C Concordance
> "Take all coaccess > 0.5 connections; overlap them with HiCCUPS loop calls from a published Hi-C dataset matching the cell type. Report the % of Cicero connections supported by Hi-C."

### SCENIC+ TF Networks
> "Run SCENIC+ on Multiome data: pycisTopic for ATAC topics, motif enrichment per topic (cisTarget), then SCENIC+ to integrate with paired RNA. Output TF-enhancer-gene triplets."

## What the Agent Will Do

1. Verify input is scATAC peak matrix (from Signac/ArchR/SnapATAC2)
2. Reduce dimensions (LSI + UMAP) -- typically already done
3. Build metacell-aggregated CDS (Cicero) or use ArchR's k-NN aggregation
4. Run co-accessibility within `window=500000` (cis) by default
5. Apply graphical lasso with default or estimated alpha
6. Filter to high-confidence connections (coaccess > 0.25 standard; > 0.5 stringent)
7. Optionally per-cluster runs to capture cell-type specificity
8. Map connections to genes via promoter overlap (TSS +/- 2 kb)
9. Compare against Hi-C / Micro-C if reference available
10. Output the strong-pair table (Peak1, Peak2, coaccess) for visualization

## Tool Decision Quick Reference

| Setting | Tool |
|---------|------|
| ATAC-only standard | Cicero |
| ATAC-only inside ArchR | ArchR getCoAccessibility |
| Multiome (RNA + ATAC) direct linking | LinkPeaks (Signac) |
| Multiome TF networks | SCENIC+ |
| Pre-computed reference | GeneHancer / FANTOM5 / EpiMap |

## Connection Score Quick Reference

| Score | Use |
|-------|-----|
| > 0.5 | Stringent / high-confidence |
| > 0.25 | Standard reporting (Cicero vignette convention) |
| > 0.05 | Exploratory; many false positives |
| 0.05 or below | Not biologically meaningful |
| < 0 | Anti-co-variation, dropped by every positive cutoff |
| NA | Unscored pair; dropped before counting |

## Tips

- Co-accessibility is NOT 3D contact. It's a statistical association from cell-to-cell co-variation. Hi-C concordance is partial and dataset dependent.
- Run Cicero per-cluster for heterogeneous datasets; pooled metacells dilute cell-type-specific connections.
- Default `window=500000` (500 kb cis) excludes most distal regulation. Widen `window` for gene desertless TADs; trans connections require Hi-C, not co-accessibility.
- Cicero's alpha parameter dramatically shifts results; use `estimate_distance_parameter()` for data-driven choice.
- Multiome data enables direct enhancer-gene correlation (LinkPeaks) which is more direct than co-accessibility alone. Restrict `genes.use` to a gene set: 150 chr1 genes took 8 min on 2,487 cells, all 1,530 chr1 genes did not finish in 25 min.
- SCENIC+ is multi-step and complex; budget 1-2 days for setup. Outputs are TF -> enhancer -> gene triples.
- Cicero output lists every pair in both orientations; the script keeps one row per unordered pair. It models only the chromosomes present in the peak set (chr1 1-30 Mb, 1,726 peaks x 3,277 cells: 8.5-14 min on 24 cores).
- Enhancer-gene mapping needs a TSS BED (chr, start, end, gene name in column 4; e.g. GENCODE protein-coding TSS) passed as the 4th CLI argument or `tss_bed`. Rows with `enhancer_is_promoter=TRUE` are promoter-promoter pairs.
- Co-accessibility predictions are hypotheses. Validate with Hi-C, ChIP-seq, or CRISPRi-FlowFISH.
- Heterochromatic 3D contacts are constitutive (don't vary in accessibility) and won't appear in co-accessibility. This is expected, not a bug.
- For published reference enhancer-gene pairs: GeneHancer (cell-type-agnostic), FANTOM5 (CAGE-based), EpiMap (~833 epigenomes).
- The same biology gives different connection counts at different metacell k values; lower k = more variability captured but slower runtime.

## Related Skills

- atac-seq/single-cell-atac - scATAC preprocessing (input)
- atac-seq/consensus-peakset - Peak set for connection inference
- atac-seq/motif-deviation - chromVAR (complement)
- gene-regulatory-networks/scenic-regulons - Standalone SCENIC for TF networks
- hi-c-analysis/loop-calling - Hi-C physical contacts (validation)
- hi-c-analysis/contact-pairs - Hi-C / Micro-C contact pairs
- single-cell/multimodal-integration - Multiome integration
- chip-seq/peak-annotation - Cross-validate with ChIP
- pathway-analysis/gsea - Downstream gene-level enrichment
