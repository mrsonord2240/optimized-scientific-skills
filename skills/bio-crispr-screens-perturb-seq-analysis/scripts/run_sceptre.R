# Purpose: SCEPTRE low-MOI differential expression (NB GLM + conditional resampling) end to end.
# Input:   an .rds holding a list with response_matrix, grna_matrix, grna_target_data_frame, extra_covariates
#          (data frame, may be NULL) and discovery_pairs (data frame); or `--example` to run sceptre's own
#          bundled lowmoi_example_data with construct_trans_pairs().
# Output:  TSV of per-gene-per-perturbation results (p_value, log_2_fold_change, significant, ...) at out_tsv.
# Usage:   Rscript scripts/run_sceptre.R input.rds results.tsv
#          Rscript scripts/run_sceptre.R --example results.tsv
args <- commandArgs(trailingOnly = TRUE)
in_arg <- args[1]
out_tsv <- args[2]
library(sceptre)

if (in_arg == '--example') {
  data(lowmoi_example_data, package = 'sceptre')
  inp <- lowmoi_example_data
  inp$discovery_pairs <- NULL
} else {
  inp <- readRDS(in_arg)
}
response_matrix <- inp$response_matrix
grna_matrix <- inp$grna_matrix
grna_target_data_frame <- inp$grna_target_data_frame
covariates_df <- inp$extra_covariates
pairs_df <- inp$discovery_pairs

# Input: sce object or sparse matrix + metadata
# Required: gene_expression_matrix, perturbation_indicator (binary per cell per pert),
#           technical_factors (batch, n_genes, etc.)

# For each gene + perturbation pair:
# Current sceptre API is a pipeline of composable steps:
sceptre_object <- import_data(response_matrix, grna_matrix, grna_target_data_frame,
                              moi = 'low', extra_covariates = covariates_df)
if (is.null(pairs_df)) pairs_df <- construct_trans_pairs(sceptre_object)
sceptre_object <- set_analysis_parameters(sceptre_object, discovery_pairs = pairs_df)
sceptre_object <- assign_grnas(sceptre_object)
sceptre_object <- run_qc(sceptre_object)
sceptre_object <- run_calibration_check(sceptre_object)
sceptre_object <- run_discovery_analysis(sceptre_object)
results <- get_result(sceptre_object, analysis = 'run_discovery_analysis')
# Output: per-gene-per-pert p-value, log-fold-change, FDR
write.table(results, out_tsv, sep = '\t', quote = FALSE, row.names = FALSE)
