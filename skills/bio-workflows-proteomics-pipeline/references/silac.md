# SILAC workflow

Read this when the data are SILAC ratios from MaxQuant `proteinGroups.txt` (moderated one-sample limma on log2 H/L ratios).

### SILAC Workflow
Caveat: heavy-Arg -> heavy-Pro metabolic conversion biases ratios for proline-containing peptides (under-counts the heavy channel), and labeling efficiency must be checked (residual light reads as down-regulation). Route to proteomics/quantification for the mechanics.
Run `scripts/silac_limma.R proteinGroups.txt silac_results.csv` (non-finite ratios to NA, >= 2 finite ratios per protein, intercept-only moderated limma, BH-adjusted p).
