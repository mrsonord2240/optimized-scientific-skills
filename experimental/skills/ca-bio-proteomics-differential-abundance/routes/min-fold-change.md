# Minimum fold change

```bash
Rscript scripts/limma_de.R proteinGroups.txt samples.csv Treatment-Control results.csv normalize=median fold=1.2
```

- `fold=1.2` is 1.2-fold; use the threshold the request gives (1.5-fold, 2-fold). The script runs `treat()` and writes `topTreat` columns: `adj.P.Val` tests `|log2FC| > log2(fold)`.
- Do not use `topTable(lfc=)` or a volcano double filter in its place.
- Input and other options are those of `routes/protein-table.md`.

Done when `results.csv` exists and you state the fold-change threshold with the count.
