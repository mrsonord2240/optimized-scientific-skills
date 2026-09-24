## SCEPTRE for Low-MOI Differential Expression

**Why this matters:** Standard differential-expression tools (DESeq2, MAST) assume Gaussian-mixture distribution and fail at single-cell scale with sparse, zero-inflated data. SCEPTRE (Katsevich Lab, 2021; low-MOI variant Barry 2024 Genome Biol) uses a negative-binomial GLM with conditional resampling:

1. Per gene, fit NB GLM: `log(expr_g) ~ pert_indicator + technical_factors`
2. Compute z-score for the perturbation coefficient
3. Resample the pert_indicator (conditional on counts) 500-1000 times; compute permutation null
4. Get FDR via permutation; not parametric

Pipeline (current sceptre API, composable steps): `import_data` -> `set_analysis_parameters` -> `assign_grnas` ->
`run_qc` -> `run_calibration_check` -> `run_discovery_analysis` -> `get_result`. Inputs: response (gene x cell) matrix,
gRNA matrix, `grna_target_data_frame`, technical covariates (batch, n_genes, etc.), `discovery_pairs`.

```bash
Rscript scripts/run_sceptre.R input.rds results.tsv     # input.rds: list(response_matrix, grna_matrix, grna_target_data_frame, extra_covariates, discovery_pairs)
Rscript scripts/run_sceptre.R --example results.tsv     # sceptre's bundled lowmoi_example_data
```
Output: per-gene-per-perturbation table (`p_value`, `log_2_fold_change`, `significant`, ...).

**Advantage over MAST:** SCEPTRE's permutation NB GLM is the only method that maintains calibrated FDR in pooled-screen scRNA-seq (Barry 2024 benchmark). MAST and Wilcoxon are over-confident due to data sparsity.
