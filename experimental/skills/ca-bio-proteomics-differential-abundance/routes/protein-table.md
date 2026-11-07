# Protein table, limma

`samples.csv`: columns `sample`, `condition`, optionally `batch`. The contrast is `Test-Reference`, two values of `condition`.

Raw search-engine table (MaxQuant proteinGroups-style TSV with linear `Intensity <sample>` columns):

```bash
Rscript scripts/limma_de.R proteinGroups.txt samples.csv Treatment-Control results.csv normalize=median
```

Log2 matrix CSV (first column protein ID, NA = missing), already normalized:

```bash
Rscript scripts/limma_de.R matrix.csv samples.csv Treatment-Control results.csv
```

- The script drops `Reverse`, `Potential contaminant` and `Only identified by site` rows, takes log2, filters, fits `eBayes(trend = TRUE, robust = TRUE)` and writes `results.csv` (`adj.P.Val` is the BH p) and `results_dropped.csv`.
- `<sample>` in `samples.csv` must match the text after the column prefix. Other prefix, e.g. LFQ: `prefix="LFQ intensity "`. LFQ and directLFQ columns are already normalized: leave `normalize` at its default.
- `normalize=median` for raw `Intensity` columns. `normalize=none` only when a global shift is the expected biology, and say so.
- `min_valid=2` is the observed values required in every group; raise it (3 of 4 is common) when asked.
- Paired or blocked design: put the subject in `batch`. Proteins whose values fall in different blocks per condition are dropped as non-estimable and listed.
- Report: samples per group, proteins removed as flagged, dropped by the valid-value filter, dropped as non-estimable, tested, significant.

Done when `results.csv` and `results_dropped.csv` exist and you have answered from them.

To keep the fit for other calculations: `source('scripts/limma_de.R')`, then `out <- run_limma_de(matrix, samples, 'Treatment-Control')`; `out$fit2` is the eBayes fit.
