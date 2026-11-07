# Volcano from a result table

```bash
Rscript scripts/volcano_csv.R results.csv volcano.pdf
```

CSV requires unique nonempty `gene`, `log2FoldChange`, and `padj` columns. Supply appropriate effect estimates, preferably a documented shrunken estimate for DESeq2 low-count results. The source function uses FDR 0.05, absolute LFC threshold 1, categorical Up/Down/NS, top ten significant labels ranked by `-log10(padj) * abs(LFC)`, and no y cap. Optional third argument gives a positive y cap; capped points are triangles. Missing adjusted p-values are excluded and counted; the wrapper rejects zero, negative or nonfinite nonmissing p-values and nonfinite effects. Output can overwrite an existing file. A horizontal line at `-log10(fdr)` is the adjusted-p threshold; both significance and LFC conditions determine color. Inspect labels, unplotted counts and capped points.

Done when input provenance, thresholds, missingness, label selection, cap and output are verified.
