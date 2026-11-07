# Consensus peak comparison

```bash
Rscript scripts/diff_accessibility.R samples.csv results lib 0.05 1 --case=treated --control=control --txdb=none
```

The CSV requires SampleID, Condition, Replicate, bamReads, Peaks and PeakCaller; prefilter peak files for blacklists. The command uses DESeq2, full-library normalization, min-overlap 2, summit +/- 250, FDR 0.05 and log2 fold threshold 1. Change `--method=edger` for the edgeR backend; `--txdb=none` writes unannotated CSV, while the default is hg38 TxDb. Use supported DiffBind design columns for covariates; see `references/design.md`. Positive fold means case over control. Source script's `fold` parameter tests an effect-size threshold; do not describe it as merely post-filtering.

Done when opened and closed BEDs and diagnostics are checked, and contrast, counting, normalization, thresholds and annotation are reported.
