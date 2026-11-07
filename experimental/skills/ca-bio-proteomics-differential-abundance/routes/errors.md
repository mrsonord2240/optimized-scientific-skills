# Errors and results that look wrong

| Symptom | Cause | Fix |
|---------|-------|-----|
| `median log2FC = ...: the matrix does not look normalized` | Raw intensities, log2 but never normalized | Rerun with `normalize=median`; `normalize=none` only for an expected global shift |
| `these look like linear intensities, not log2` | A CSV of linear values given as a log2 matrix | Give the raw `.tsv`, or take log2 first |
| `Columns not in the table` | Prefix plus sample name matches no column | Set `prefix=` (e.g. `prefix="LFQ intensity "`) and match `samples.csv` names |
| `prior.weights contain NA values` or `Partial NA coefficients` | Rows with no valid value, or values in different blocks per condition | `limma_de.R` filters both; read `results_dropped.csv` |
| `makeContrasts: object 'X' not found` | Contrast names are not values of `condition` | Write `Test-Reference` with exact condition values |
| Many calls, all in one direction, with `normalize=none` | Loading difference tested as biology | `normalize=median` |
| DEqMS stops on counts, or `sca.*` columns missing | Missing or zero count; limma columns read | Supply a count of at least 1 per protein; read `sca.adj.pval` |
| proDA: `object 'conditionControl' not found` | Coefficients are `Intercept`, `conditionTreatment` | The script tests the right one; in your own code test the non-reference level |
| msqrob2: `adjPval = NA` rows | No observation in a condition, or not fittable | Report from `results_undetected.csv` |
| msqrob2: ridge refused | `ridge = TRUE` on a two-group design | `ridge = FALSE` |
| msqrob2: `Variable sample is not found` | `sample` missing from `colData`, `feature` from `rowData` | Use `scripts/msqrob2_de.R` |
| `readQFeatures: 'colData' must contain a column called 'quantCols'` | QFeatures 1.16 | Add `quantCols = runs` to `colData` |
| `pe[['peptideLog']]$condition` is NULL | `colData` is on the QFeatures object | `colData(pe)$condition` |
| MSstats rows with infinite `log2FC` | `issue == oneConditionMissing` | Report as undetected, never as a ratio |
| Flag filter on `evidence.txt` drops every row | Empty flag cells are NA | Filter with `!(ev$Reverse %in% '+')` |
| Many low-\|logFC\| calls at feature level, median logFC far from 0 | Per-run median normalization of the peptide table | `references/feature_level.md` |
| Volcano with rigid streaks of on/off proteins | Imputation upstream | Rerun on the unimputed table with `routes/heavy-missingness.md` |
| Python: no protein has 2 values per group | A group has one sample | Not a two-sample test; do not run it |
