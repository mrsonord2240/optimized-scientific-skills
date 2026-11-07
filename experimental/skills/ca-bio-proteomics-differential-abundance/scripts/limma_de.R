#!/usr/bin/env Rscript
# Purpose: limma differential abundance from a protein table: drop flagged rows, log2, valid-value filter,
#   design with batch as a covariate, estimability filter, empirical-Bayes moderation (trend + robust), BH
#   table, and a centring check that stops when the matrix was not normalized.
# Inputs:  protein table: raw search-engine TSV (MaxQuant proteinGroups-style, linear `Intensity <sample>`
#          columns) or a log2 matrix CSV (first column = protein ID, NA = missing);
#          samples CSV (columns: sample, condition, optionally batch); contrast as levels of `condition`
#          written Treatment-Control; output CSV.
# Usage:   Rscript limma_de.R proteins.tsv samples.csv Treatment-Control results.csv [key=value ...]
#            normalize=check|median|none  check (default) stops if the median log2FC shows an un-normalized
#                                matrix; median centres every sample first; none fits as given
#            min_valid=2         observed values needed in every group
#            prefix="Intensity " raw TSV only: column prefix before the sample name (e.g. "LFQ intensity ")
#            id=<column>         raw TSV only: protein ID column (default Protein IDs)
#            format=raw|log2     default: .tsv/.txt = raw, .csv = log2
#            fold=1.2            test against a minimum fold change with treat() instead of against zero
#            shrink=ashr         add logFC_shrunk and lfsr columns
#          or source() this file and call run_limma_de() to keep `fit2`.
# Output:  results CSV (protein, logFC, AveExpr, t, P.Value, adj.P.Val [, B]) plus <results>_dropped.csv
#          listing every protein not tested and why.
# Checked: limma 3.62.2 (R 4.4.3 / Bioconductor 3.20).
suppressPackageStartupMessages(library(limma))
.here <- local({
  f <- grep('^--file=', commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub('^--file=', '', f[1]))) else getwd()
})
source(file.path(.here, 'protein_input.R'))

OFFSET_MAX <- 0.1  # |median log2FC| trip-wire for a protein matrix; empirical, not a distributional bound

run_limma_de <- function(protein_matrix, sample_info, contrast, min_valid = 2,
                         normalize = c('check', 'median', 'none'), fold = NULL) {
  normalize <- match.arg(normalize)
  cond <- factor(sample_info$condition)  # contrast names must be its levels
  sample_info$condition <- cond
  protein_matrix <- protein_matrix[, sample_info$sample, drop = FALSE]
  if (normalize == 'median') protein_matrix <- center_medians(protein_matrix)

  # Valid-value filter BEFORE lmFit: >= min_valid values in every group (study choice; 3 of 4 is common)
  n_valid <- sapply(levels(cond), function(g) rowSums(!is.na(protein_matrix[, cond == g, drop = FALSE])))
  n_valid <- matrix(n_valid, ncol = nlevels(cond), dimnames = list(rownames(protein_matrix), levels(cond)))
  keep_valid <- apply(n_valid >= min_valid, 1, all)
  dropped_valid <- data.frame(protein = rownames(protein_matrix)[!keep_valid], reason = rep('too few observed values', sum(!keep_valid)),
                              n_valid[!keep_valid, , drop = FALSE], check.names = FALSE, row.names = NULL)
  protein_matrix <- protein_matrix[keep_valid, , drop = FALSE]

  if ('batch' %in% names(sample_info)) {  # batch in the model, not removed first
    sample_info$batch <- factor(sample_info$batch)
    design <- model.matrix(~0 + condition + batch, data = sample_info)
  } else {
    design <- model.matrix(~0 + condition, data = sample_info)
  }
  colnames(design)[seq_len(nlevels(cond))] <- levels(cond)

  fit <- lmFit(protein_matrix, design)
  # Estimability filter: the per-group count does not make the contrast estimable under a blocked design
  # ('Partial NA coefficients for N probe(s)'). Keep only fully estimated rows with residual df.
  estimable <- fit$df.residual > 0 & rowSums(is.na(fit$coefficients)) == 0
  dropped_nonestimable <- data.frame(protein = rownames(fit)[!estimable],
                                     reason = rep('not estimable in this design', sum(!estimable)))  # rep(): a zero-row frame cannot recycle a scalar
  fit <- fit[estimable, ]

  contrast_matrix <- makeContrasts(contrasts = contrast, levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)  # trend mandatory for label-free; robust Winsorizes outliers
  results <- topTable(fit2, coef = 1, number = Inf, adjust.method = 'BH')  # adj.P.Val is the BH p

  # Centring check: most proteins are unchanged, so the median log2FC of a normalized matrix is near 0.
  # A loading difference between the groups shifts every protein and turns into calls in one direction.
  offset <- median(results$logFC, na.rm = TRUE)
  if (normalize == 'check' && abs(offset) > OFFSET_MAX) stop(sprintf(paste(
    'median log2FC = %+.3f: the matrix does not look normalized, so a per-sample loading difference is',
    'being tested as biology. Rerun with normalize=median, or with normalize=none if a global shift is',
    'the expected biology (then say so in the report).'), offset))

  fit_treat <- NULL
  if (!is.null(fold)) {  # treat() re-estimates the prior; its trend/robust default to FALSE
    fit_treat <- treat(fit2, lfc = log2(fold), trend = TRUE, robust = TRUE)
    results <- topTreat(fit_treat, coef = 1, number = Inf)  # no B column
  }
  list(fit2 = fit2, fit_treat = fit_treat, results = results, design = design, offset = offset,
       dropped = rbind(dropped_valid[, c('protein', 'reason')], dropped_nonestimable), dropped_valid = dropped_valid)
}

if (sys.nframe() == 0) {  # run as a script, not source()d
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 4) stop('Usage: Rscript limma_de.R proteins.tsv samples.csv Treatment-Control results.csv [key=value ...]')
  kv <- parse_kv(args[-(1:4)])
  samples <- read_samples(args[2])
  inp <- read_protein_matrix(args[1], samples, kv)
  normalize <- if (!is.null(kv$normalize)) kv$normalize else 'check'
  fold <- if (!is.null(kv$fold)) as.numeric(kv$fold) else NULL
  out <- run_limma_de(inp$matrix, samples, args[3], if (!is.null(kv$min_valid)) as.integer(kv$min_valid) else 2L,
                      normalize, fold)
  res <- out$results
  if (identical(kv$shrink, 'ashr')) {
    if (!requireNamespace('ashr', quietly = TRUE)) stop('shrink=ashr needs the ashr package')
    fit2 <- out$fit2
    ok <- !is.na(fit2$coefficients[, 1]) & !is.na(fit2$s2.post)  # ash() gives PosteriorMean 0 for NA rows
    se <- sqrt(fit2$s2.post[ok]) * fit2$stdev.unscaled[ok, 1]
    sh <- ashr::ash(fit2$coefficients[ok, 1], se, mixcompdist = 'normal')$result
    idx <- match(rownames(res), rownames(fit2$coefficients)[ok])
    res$logFC_shrunk <- sh$PosteriorMean[idx]
    res$lfsr <- sh$lfsr[idx]
  }
  write.csv(cbind(protein = rownames(res), res), args[4], row.names = FALSE)
  dropped <- rbind(inp$removed, out$dropped)
  write.csv(dropped, sub('\\.csv$', '_dropped.csv', args[4]), row.names = FALSE)
  cat('normalize:', normalize, '| median log2FC:', sprintf('%+.3f', out$offset), '\n')
  cat('flagged rows removed:', nrow(inp$removed), '| dropped by valid-value filter:', nrow(out$dropped_valid),
      '| dropped as non-estimable:', nrow(out$dropped) - nrow(out$dropped_valid), '| tested:', nrow(res), '\n')
  cat('significant (adj.P.Val < 0.05)', if (!is.null(fold)) sprintf('against a %.2f-fold minimum', fold) else '',
      ':', sum(res$adj.P.Val < 0.05), '\n')
}
