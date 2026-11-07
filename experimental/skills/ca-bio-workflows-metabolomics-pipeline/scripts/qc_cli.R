args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop('usage: Rscript scripts/qc_cli.R input.rds output.rds')
if (file.exists(args[2])) stop('output already exists')
x <- readRDS(args[1])
feat <- x$feat; defs <- x$defs
if (!is.null(x$pd)) {
  stopifnot(all(c('sample_group', 'injection_order', 'batch') %in% names(x$pd)))
  sample_class <- x$pd$sample_group; injection_order <- x$pd$injection_order; batch_id <- x$pd$batch
} else {
  sample_class <- x$sample_class; injection_order <- x$injection_order; batch_id <- x$batch_id
}
source('scripts/stage2_qc_impute.R')
saveRDS(list(imputed=imputed, sample_class=sample_class, defs=defs, normalized=normalized), args[2])
