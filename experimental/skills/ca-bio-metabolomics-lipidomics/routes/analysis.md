# Differential analysis and enrichment

```bash
Rscript examples/lipidomics_workflow.R
```

This source example uses lipidr's bundled normalized dataset and a fixed contrast, writing to a temporary directory. It demonstrates differential analysis, a volcano and class/chain enrichment; it is not a user-data CLI. Adapt import, sample annotation, normalization, contrast and output paths for the actual cohort before use. Inspect unclassified names and report tested sets, correction, effect sizes and resolution levels. The source's fixed thresholds are adjusted p < 0.05 and absolute log2 fold change > 1. See [interpretation](../references/interpretation.md).

Done when input provenance, contrast, normalization, uncertainty and results are verified for the actual study.
