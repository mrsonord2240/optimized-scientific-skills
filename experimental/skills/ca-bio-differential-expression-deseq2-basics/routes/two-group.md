# Two groups

```bash
Rscript scripts/deseq2_de.R counts.csv samples.csv results.csv ref=control
```

`counts.csv`: first column gene ID, one column per sample, raw integer counts. `samples.csv`: first column sample ID, then a `condition` column (another name: `condition=<column>`).

- `ref=` is the denominator. `log2FoldChange` is the other level over `ref`.
- Significant genes are the rows with `padj < 0.05` in `results.csv`.
- Fewer than 3 samples in a group: run it, and report the result as exploratory.
- Report the counts the script prints: samples per level, genes tested, significant, and `padj = NA` by cause.

Done when `results.csv` exists and you have answered from it.

Only if it applies:

- More than half the genes change in one direction (bacterial stress, viral shutoff, global amplification): `references/size-factors.md`.
- Bacterial or archaeal data: `references/prokaryotic-rna-seq.md`.
