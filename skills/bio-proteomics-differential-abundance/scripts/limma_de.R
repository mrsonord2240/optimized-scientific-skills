#!/usr/bin/env Rscript
# Purpose: limma differential abundance on a log2, normalized protein matrix: valid-value filter, design with
#   batch as a covariate, estimability filter, empirical-Bayes moderation (trend + robust), BH table.
# Inputs:  matrix CSV (first column = protein ID, one column per sample, log2 intensities, NA = missing);
#          sample CSV (columns: sample, condition, and optionally batch/donor named 'batch');
#          contrast as levels of `condition`, e.g. Treatment-Control; output CSV.
# Usage:   Rscript limma_de.R matrix.csv samples.csv Treatment-Control results.csv [min_valid=2]
#          or source() this file and call run_limma_de() to keep `fit2` for treat(), DEqMS and ashr.
# Checked: limma 3.62.2 (R 4.4.3 / Bioconductor 3.20).
suppressPackageStartupMessages(library(limma))

run_limma_de <- function(protein_matrix, sample_info, contrast, min_valid = 2) {
  cond <- factor(sample_info$condition)  # first level is the reference; contrast names must be its levels
  sample_info$condition <- cond
  protein_matrix <- protein_matrix[, sample_info$sample, drop = FALSE]

  # Valid-value filter BEFORE lmFit: >= min_valid values in every group (study choice; 3 of 4 is common)
  n_valid <- sapply(levels(cond), function(g) rowSums(!is.na(protein_matrix[, cond == g, drop = FALSE])))
  keep_valid <- apply(n_valid >= min_valid, 1, all)
  dropped_valid <- rownames(protein_matrix)[!keep_valid]  # report these (e.g. undetected in one group)
  protein_matrix <- protein_matrix[keep_valid, , drop = FALSE]

  if ('batch' %in% names(sample_info)) {  # batch in the model, not removed first
    sample_info$batch <- factor(sample_info$batch)
    design <- model.matrix(~0 + condition + batch, data = sample_info)
  } else {
    design <- model.matrix(~0 + condition, data = sample_info)
  }
  colnames(design)[seq_len(nlevels(cond))] <- levels(cond)

  fit <- lmFit(protein_matrix, design)
  # Estimability filter: the per-group count above does not make the contrast estimable under a blocked
  # design ('Partial NA coefficients for N probe(s)'). Keep only fully estimated rows with residual df.
  estimable <- fit$df.residual > 0 & rowSums(is.na(fit$coefficients)) == 0
  dropped_nonestimable <- rownames(fit)[!estimable]  # report these; non-estimable, not "not significant"
  fit <- fit[estimable, ]

  contrast_matrix <- makeContrasts(contrasts = contrast, levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)  # trend mandatory for label-free; robust Winsorizes outliers

  results <- topTable(fit2, coef = 1, number = Inf, adjust.method = 'BH')
  # columns: logFC, AveExpr, t, P.Value, adj.P.Val, B  (adj.P.Val is the BH p; there is no $FDR)
  list(fit2 = fit2, results = results, design = design,
       dropped_valid = dropped_valid, dropped_nonestimable = dropped_nonestimable)
}

if (sys.nframe() == 0) {  # run as a script, not source()d
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 4) stop('Usage: Rscript limma_de.R matrix.csv samples.csv Treatment-Control results.csv [min_valid=2]')
  m <- as.matrix(read.csv(args[1], row.names = 1, check.names = FALSE))
  s <- read.csv(args[2], stringsAsFactors = FALSE)
  min_valid <- if (length(args) >= 5) as.integer(args[5]) else 2L
  out <- run_limma_de(m, s, args[3], min_valid)
  write.csv(out$results, args[4])
  cat('tested:', nrow(out$results), '| dropped by valid-value filter:', length(out$dropped_valid),
      '| dropped as non-estimable:', length(out$dropped_nonestimable), '\n')
  cat('significant (adj.P.Val < 0.05):', sum(out$results$adj.P.Val < 0.05), '\n')
}
