---
name: bio-atac-seq-differential-accessibility
description: Identify differentially accessible chromatin regions across conditions using DiffBind, csaw, DESeq2, or edgeR. Use when comparing ATAC-seq accessibility between treatment groups, choosing between consensus-peak vs sliding-window approaches, picking the correct normalization (full library vs reads-in-peaks), correcting batch with SVA/RUVseq, or interpreting log2FC and FDR thresholds in a chromatin context.
tool_type: r
primary_tool: DiffBind
license: MIT
category: Data Analysis
author: GPTomics
---

# Differential Accessibility

Find chromatin regions that change accessibility between conditions. Build a sample-by-region count matrix, normalize for library size and chromatin compaction, fit a generalized linear model (negative-binomial), and extract regions with significant accessibility change.

## Workflow

1. Prepare input: per-sample BAM files (deduplicated, chrM-stripped, MAPQ-filtered), per-sample peak files (narrowPeak format), and a sample sheet defining conditions and any covariates (batch, donor, time).
2. Choose the analysis strategy:
   - Consensus-peak + DiffBind (default, ATAC-aware): stable peaks across conditions, 2-3+ reps per group
   - csaw sliding windows: peak structure differs dramatically between conditions, peak-free analysis
   - Direct DESeq2/edgeR on existing peak-count matrix: for existing featureCounts output
3. Select normalization approach: full-library (script default, preserves global shifts), reads-in-peaks or `native` RLE/TMM (when background varies independently of biology), or spike-in (when global compaction is biological; recipe only).
4. Add batch/donor/time covariates to the design when applicable.
5. Apply SVA or RUVseq for hidden batch effects when empirical surrogates are needed. The script's `--sva` mode fits DESeq2 directly: it ignores `--method` and `--design`, skips the blacklist filter and heatmap, and applies `lfc_thr` as a hard post-filter, whereas the default path passes it to DiffBind as an lfcThreshold test (stricter), so hit counts are not comparable across modes. With few samples `n.sv` is capped (4 samples allow 1) and surrogates can absorb real signal.
6. Report results at FDR threshold and effect-size cutoff; annotate differential peaks to genes.
7. Document the consensus peakset strategy, normalization choice, tool versions, and software environment.

## Routed material

- Run DiffBind (DESeq2 or edgeR backend, optional covariate design or surrogate variables) with [`scripts/diff_accessibility.R`](scripts/diff_accessibility.R); its header documents every CLI option. Verify sample-sheet format, peak files, BAM paths and Condition labels first. Only this route is shipped as executable code; edgeR-on-a-matrix, csaw, RUVseq and spike-in are recipes in the references.
- Use [`references/method-reference.md`](references/method-reference.md) for algorithmic comparison, decision trees by scenario, normalization strategy, failure modes, reconciliation patterns, and references.
- Use [`references/usage-guide.md`](references/usage-guide.md) for quick-start examples, method selection, and tool tips.
- Refer to related Skills for consensus peakset generation, QC, single-cell ATAC, or co-accessibility analysis (see Related Skills below).

## Inputs and outputs

Inputs: Per-replicate deduplicated BAM files with chrM-stripped reads; per-replicate peak files (narrowPeak); sample sheet CSV with SampleID, Condition, Replicate, bamReads, Peaks, PeakCaller, and optional Tissue/Factor/Treatment covariate columns (the only extra names DiffBind designs accept).

Outputs: Differentially accessible peak BED files (opened.bed, closed.bed); annotated CSV with peak coordinates and nearest genes; diagnostic plots (PCA, MA, volcano, heatmap, annotation distribution).

## Version compatibility

Reference examples tested with: DiffBind 3.12+, DESeq2 1.42+, edgeR 4.0+, csaw 1.36+, limma 3.58+, GenomicRanges 1.54+, ChIPseeker 1.38+, sva 3.50+, RUVSeq 1.36+, Subread 2.0+ (featureCounts). Before using code patterns, verify installed versions with `packageVersion('<pkg>')` in R or `conda list` in conda environments. Inspect package help (`?function_name`) to confirm parameter availability.

## License and provenance

This derived Skill retains GPTomics/bioSkills material by Domen Jemec, from commit `d91ed3d563019e649dc854c56ccd62551359488a`. See [`LICENSE`](LICENSE) for the preserved MIT notice.

## Related Skills

- atac-seq/atac-peak-calling - Generate per-replicate peaks
- atac-seq/consensus-peakset - Build differential-ready consensus peakset
- atac-seq/atac-qc - Pre-screen and drop failing replicates
- atac-seq/single-cell-atac - Pseudobulk-level differential per cluster
- atac-seq/co-accessibility - Identify cis-regulatory connections among DA peaks
- differential-expression/deseq2-basics - Underlying DESeq2 patterns
- differential-expression/de-results - Effect-size reporting and shrinkage
- chip-seq/differential-binding - Same DiffBind workflow, ChIP context
- pathway-analysis/go-enrichment - Downstream gene-level enrichment of DA-associated genes
