#!/usr/bin/env Rscript
# SNP-by-trait beta and SE TSVs: first column SNP ID, remaining columns trait names.
# Usage: Rscript scripts/hyprcoloc_traits.R betas.tsv ses.tsv results.tsv
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop('usage: hyprcoloc_traits.R betas.tsv ses.tsv results.tsv')
suppressPackageStartupMessages(library(hyprcoloc))
read_matrix <- function(path) {
  x <- read.delim(path, row.names = 1, check.names = FALSE)
  if (!nrow(x) || !ncol(x) || anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)))
    stop('each matrix needs unique SNP and trait names and at least one row and column')
  if (!all(vapply(x, is.numeric, logical(1)))) stop('matrix values must be numeric')
  x <- as.matrix(x)
  if (any(!is.finite(x))) stop('matrix contains missing or non-finite values')
  x
}
betas <- read_matrix(args[1])
ses <- read_matrix(args[2])
if (!identical(dimnames(betas), dimnames(ses))) stop('SNP and trait names and order must match')
if (any(ses <= 0)) stop('standard errors must be positive')
res <- hyprcoloc(effect.est = betas, effect.se = ses, trait.names = colnames(betas),
                  snp.id = rownames(betas), reg.thresh = 0.7, align.thresh = 0.7)
if (is.null(res$results)) stop('HyPrColoc returned no results table; inspect installed API')
write.table(res$results, args[3], sep = '\t', quote = FALSE, row.names = FALSE)
