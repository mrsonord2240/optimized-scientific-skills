# Values for PCA, heatmaps or machine learning

```r
source('scripts/deseq2_de.R')
out <- run_deseq2(dds, condition = 'condition', ref = 'control')   # or reuse out from your run
vsd <- vst(out$dds, blind = FALSE)
plotPCA(vsd, intgroup = 'condition')
write.csv(assay(vsd), 'vst_values.csv')
```

- Use these values for plots, clustering and model features only. Never feed them back into a differential expression test.
- Fewer than 1,000 genes: `varianceStabilizingTransformation(out$dds, blind = FALSE)`.
- Fewer than 30 samples whose library sizes differ more than four-fold: `rlog(out$dds, blind = FALSE)`.

Done when the plot or the values file is written.

More detail: `references/transforms.md`.
