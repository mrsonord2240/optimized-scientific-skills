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
# Recommended: runs R in a child process, validates the TSV, and writes results.tsv.run.json.
# Pass --rscript with a pinned project R wrapper when your environment needs one.
python scripts/run_sceptre_safe.py input.rds results.tsv
python scripts/run_sceptre_safe.py --example results.tsv
```
Output: per-gene-per-perturbation table (`p_value`, `log_2_fold_change`, `significant`, ...).

The safe wrapper refuses a pre-existing output path, and returns success only if the worker exits successfully and produces a nonempty structurally valid table. The manifest records worker status and validation; a nonzero worker exit is a failure even if the worker wrote a table. Pin and retest the R/sceptre stack before production analysis. The R script is the worker, not the recommended public entrypoint.

**Advantage over MAST:** SCEPTRE's permutation NB GLM is the only method that maintains calibrated FDR in pooled-screen scRNA-seq (Barry 2024 benchmark). MAST and Wilcoxon are over-confident due to data sparsity.
