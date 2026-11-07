# DE results table

```bash
Rscript scripts/ora_go.R de=de_results.csv results.csv gene_col=gene_id keytype=SYMBOL ont=BP padj=0.05 lfc=1 direction=both
```

The table needs the gene column plus `padj`, `log2FoldChange`, `pvalue`. Foreground = `padj < padj` and `|log2FoldChange| > lfc`. Universe = rows with a non-NA `pvalue` (the genes DESeq2 actually tested), not all rows and not the genome.

- If the request states its own cutoffs, pass them as `padj=` and `lfc=`.
- `direction=up` or `direction=down` runs one direction; run both when direction matters, since mixed lists cancel.
- Other options (`ont=ALL`, `orgdb=`, `simplify=false`, `keytype=`) are those of `routes/ora.md`.
- Length-bias correction for RNA-seq: `routes/length-bias.md`.

Done when `results.csv` exists and the answer reports foreground and universe sizes from `results.csv.log`.
