# MSstats Workflow (feature-level model)

Read this when the design needs a feature-level mixed model (MSstats) instead of limma on a protein matrix: MaxQuant `evidence.txt` + `proteinGroups.txt` input, peptide/feature-level data, or a comparison matrix with more than two conditions.

## MSstats Workflow

Run `scripts/msstats_maxquant.R evidence.txt proteinGroups.txt annotation.csv msstats_comparison.csv`. Its comments carry the `quote = ''` / row-count guard, the `equalizeMedians` symmetry caveat and the contrast-matrix rules for three or more conditions.
