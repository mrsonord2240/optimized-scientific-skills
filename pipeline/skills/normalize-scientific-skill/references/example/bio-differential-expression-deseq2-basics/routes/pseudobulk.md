# Single-cell counts from several donors

```bash
Rscript scripts/pseudobulk_de.R <10x_dir> cells.tsv results.tsv ref=control
```

`<10x_dir>` holds `matrix.mtx`, `features.tsv` and `barcodes.tsv` (plain or `.gz`) with raw counts for ONE cell type. `cells.tsv`: first column barcode, then `donor` and `condition` columns (other names: `donor=<column> condition=<column>`). Several cell types: subset and run once per type.

The script sums raw counts per donor and condition, drops thin samples, fits DESeq2 with the donor as a blocking factor when donors appear in both conditions, and writes one table with gene symbols.

- The donor is the replicate, not the cell. Never test cells against cells for a condition effect.
- `min_cells=10` drops a donor-condition sample with fewer cells. If the request states its own minimum, pass it.
- `min_count=10` keeps genes with at least that total count. If the request states its own, pass it.
- If the request says to adjust p-values across all tested genes, add `all_genes_bh=true`.
- If the request defines fold change as a ratio of mean normalized expression, compute `log2(mean_norm_<test> / mean_norm_<ref>)` from the output columns. Keep `padj` from the same file.
- Read the gene list from `results.tsv` by its `symbol` column.

Done when `results.tsv` exists and you have answered from it. Do not rerun with a t-test, edgeR or limma to compare.
