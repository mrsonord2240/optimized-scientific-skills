---
name: bio-atac-seq-co-accessibility
description: Infer cis-regulatory connections (peak-to-peak co-accessibility) from scATAC-seq using Cicero, ArchR getCoAccessibility, or SCENIC+. Use when linking enhancer accessibility to promoter accessibility, identifying enhancer-gene pairs from chromatin alone (without paired RNA), running gene-regulatory inference combining ATAC + RNA, or comparing predicted regulatory contacts against Hi-C/Micro-C ground truth.
tool_type: r
primary_tool: cicero
license: MIT
category: Data Analysis
author: GPTomics
---

# Co-accessibility (cis-Regulatory Linkage)

Infer enhancer-promoter and enhancer-enhancer cis-regulatory connections from scATAC-seq using cell-to-cell variability in joint peak accessibility. Cicero, ArchR getCoAccessibility, and SCENIC+ produce peak-pair connection scores; thresholding generates enhancer-gene candidates.

**Key principle:** Co-accessibility is NOT 3D contact—it's a statistical association from cell-to-cell co-variation. Strong co-accessibility overlaps Hi-C/Micro-C contacts only partly (concordance is dataset dependent) and is not equivalent. Validate against orthogonal data (Hi-C, CRISPRi-FlowFISH).

## Quick Start

Choose your input and tool:

| Your data | Tool | Output | Status |
|-----------|------|--------|--------|
| scATAC peak matrix alone | Cicero or ArchR getCoAccessibility | Peak-pair connection scores (-1 to 1) | Cicero: shipped script, executed. ArchR: executed on a PBMC chr1 slice |
| Multiome (ATAC + RNA) | LinkPeaks (Signac) or SCENIC+ | Direct enhancer-gene links or TF-driven networks | LinkPeaks: executed on a bounded gene set, no shipped code. SCENIC+: not executed, static review only |
| Validation against Hi-C/Micro-C | Any co-accessibility + HiCCUPS loops | Overlap % as concordance measure | Snippet executed on a planted BEDPE |
| Pre-computed reference (any data) | GeneHancer, FANTOM5, EpiMap | Published enhancer-gene pairs (cell-type-agnostic) | Prose only, no code |

## Workflow

1. **Prepare input:** scATAC peak-cell matrix (from Signac, ArchR, or SnapATAC2 preprocessing) with a UMAP or dimension reduction.
2. **Select tool:** Cicero for ATAC-only, ArchR getCoAccessibility if already in ArchR, LinkPeaks for Multiome, SCENIC+ for TF networks, or reference databases for published pairs.
3. **Run co-accessibility:** Build metacells (Cicero) or use tool's aggregation, compute correlation across cis window (default 500 kb), apply graphical lasso regularization.
4. **Filter:** Threshold on connection score (> 0.25 standard, > 0.5 stringent).
5. **Map to genes (optional):** Overlap one anchor with promoters (TSS +/- 2 kb) to generate enhancer-gene candidate pairs.
6. **Validate:** Compare against Hi-C/Micro-C or experimental enhancer-promoter interaction assays.

## Routed material

- **Usage guide and tool decision trees:** See [`references/usage-guide.md`](references/usage-guide.md).
- **Detailed methodology, failure modes, and citations:** See [`references/method-reference.md`](references/method-reference.md).
- **Reference R script (Cicero):** See [`scripts/cicero_workflow.R`](scripts/cicero_workflow.R) (header documents the CLI arguments and file formats). Verify installed package versions with `packageVersion('<pkg>')` before running.

## Inputs and outputs

**Expected inputs:**
- scATAC peak-cell matrix (binary or counts; from Signac, ArchR, or SnapATAC2)
- Cell metadata (cluster, cell type, or pseudotime annotations)
- Peak metadata (row names `chr_start_end`, e.g. `chr1_10244_10510`; the script rejects other formats)
- Dimension reduction (UMAP or other) for metacell construction
- Optional: paired RNA AnnData (Multiome), published Hi-C loops for validation

**Principal outputs:**
- Peak-pair connection score DataFrame (Peak1, Peak2, coaccess score), one row per unordered pair
- High-confidence subset (thresholded by connection score)
- Enhancer-gene candidate pairs (if a TSS BED with gene names in column 4 is supplied), with promoter-promoter pairs flagged
- Concordance % with Hi-C/Micro-C (if ground truth available)

## License and provenance

This derived Skill retains GPTomics/bioSkills material from commit `d91ed3d563019e649dc854c56ccd62551359488a`. See [`LICENSE`](LICENSE) for the preserved MIT notice.

## Related Skills

- atac-seq/single-cell-atac — scATAC preprocessing (input)
- atac-seq/consensus-peakset — Peak set used for connection inference
- atac-seq/motif-deviation — chromVAR for TF activity (complement)
- atac-seq/enhancer-gene-linking — ABC, ENCODE-rE2G when Hi-C is available
- gene-regulatory-networks/scenic-regulons — Standalone SCENIC for TF networks
- hi-c-analysis/loop-calling — Physical contacts from Hi-C (validation)
- single-cell/multimodal-integration — Multiome integration
