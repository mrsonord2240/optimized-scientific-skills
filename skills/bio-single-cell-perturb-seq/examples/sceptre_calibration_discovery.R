# Reproducible low-MOI SCEPTRE smoke run: calibration results are required before discovery.
# Before running, install the pinned package into a private SCEPTRE_R_LIB as documented in SKILL.md.
sceptre_lib <- Sys.getenv("SCEPTRE_R_LIB")
if (!nzchar(sceptre_lib) || !dir.exists(sceptre_lib)) {
  stop("Set SCEPTRE_R_LIB to the private library containing sceptre before running this example.")
}
.libPaths(c(sceptre_lib, .libPaths()))
library(sceptre)

out_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out_dir) || !nzchar(out_dir)) out_dir <- "sceptre-results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data(lowmoi_example_data, package = "sceptre")
rna_counts <- lowmoi_example_data$response_matrix
grna_counts <- lowmoi_example_data$grna_matrix
grna_targets <- lowmoi_example_data$grna_target_data_frame
extra_covariates <- lowmoi_example_data$extra_covariates

obj <- import_data(response_matrix = rna_counts, grna_matrix = grna_counts,
                   grna_target_data_frame = grna_targets, moi = "low",
                   extra_covariates = extra_covariates)
pairs <- head(construct_trans_pairs(obj), 4)
obj <- set_analysis_parameters(obj, discovery_pairs = pairs)
obj <- assign_grnas(obj, method = "mixture")
obj <- run_qc(obj)
obj <- run_calibration_check(obj)
calibration <- get_result(obj, analysis = "run_calibration_check")
stopifnot(nrow(calibration) > 0)
obj <- run_discovery_analysis(obj)
results <- get_result(obj, analysis = "run_discovery_analysis")
stopifnot(nrow(results) == nrow(pairs))

write.table(calibration, file.path(out_dir, "calibration.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(results, file.path(out_dir, "discovery.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("PASS sceptre=", as.character(packageVersion("sceptre")),
    " calibration=", nrow(calibration), " discovery=", nrow(results), "\n", sep = "")
