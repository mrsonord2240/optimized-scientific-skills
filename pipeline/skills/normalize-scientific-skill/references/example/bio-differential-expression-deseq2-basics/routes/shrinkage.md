# Fold changes for ranking or plots

Add `shrink=apeglm` to the run:

```bash
Rscript scripts/deseq2_de.R counts.csv samples.csv results.csv ref=control shrink=apeglm
```

The output gains `log2FoldChange_shrunk`.

- Rank for GSEA by `log2FoldChange_shrunk`, or by `stat`. Not by the raw `log2FoldChange`, which low-count genes dominate.
- Plot the shrunken column on a volcano or MA plot.
- Significance still comes from `padj`. Shrinking changes the fold change, not the p-value.
- Do not filter on the shrunken fold change and then claim the false discovery rate covers that filter. For "changed by more than X-fold" as a tested claim, see `references/lfc-shrinkage.md`.
- If apeglm is not installed, use `shrink=ashr`.

Done when the ranked list or plot is made from the shrunken column.
