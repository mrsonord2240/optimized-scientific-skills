#!/usr/bin/env Rscript
# Purpose: the two centring checks to run before reading any feature-level (msqrob2/MSstats) result:
#   1. residual-SD ratio of the peptide table before/after per-run normalization (warning),
#   2. median log2FC of the result table (stop).
# Inputs:  peptide_wide CSV (columns 'feature', 'protein', then one intensity column per run, linear scale),
#          sample CSV (columns run, condition), normalize = TRUE|FALSE (apply per-run median centring, the
#          case the checks exist for), reference and test condition names.
# Usage:   Rscript centring_checks.R peptide_wide.csv samples.csv TRUE Control Treatment
#          or source() this file and call resid_sd(), check_sd_ratio(), check_offset() on your own objects.
# Checked: msqrob2 1.14.1, QFeatures 1.16.0, MsCoreUtils 1.18.0 (R 4.4.3 / Bioconductor 3.20).
# Trip-wires, not distributional bounds: both constants were read off ONE 4 v 4 label-free set.
OFFSET_MAX   <- 0.05  # |median log2FC| over the tested proteins
SD_RATIO_MIN <- 0.85  # residual SD after / before per-run normalization (0.76 where it produced 31 false positives)

# median within-condition SD over peptides (log2 matrix, features x runs; NA-tolerant)
resid_sd <- function(L, cond) {
  ss <- df <- 0
  for (g in levels(cond)) {
    m <- L[, cond == g, drop = FALSE]
    ss <- ss + rowSums((m - rowMeans(m, na.rm = TRUE))^2, na.rm = TRUE)
    df <- df + pmax(rowSums(!is.na(m)) - 1, 0)
  }
  median(sqrt(ss[df > 0] / df[df > 0]))
}

# 1. Residual check, only if you normalize. msqrob2 route: run this after filterFeatures and before
#    aggregateFeatures, then aggregate and test from 'peptideNorm' (aggregateFeatures(pe, i = 'peptideNorm', ...)).
#    MSstats route: fit twice, normalization = FALSE and the default, and compare median(tested$SE)
#    instead: same warning, same SD_RATIO_MIN.
check_sd_ratio <- function(before, after, cond) {
  sd_ratio <- resid_sd(after, cond) / resid_sd(before, cond)
  if (sd_ratio < SD_RATIO_MIN) warning(sprintf(
    'residual SD fell to %.2f of its un-normalized value: the normalization is removing within-condition variance and shrinking the SEs. Check the offset below and compare with normalization off.',
    sd_ratio))
  sd_ratio
}

# 2. Offset check on the result table (msqrob2: res$logFC; MSstats: tested$log2FC)
check_offset <- function(logfc) {
  offset <- median(logfc, na.rm = TRUE)
  if (abs(offset) > OFFSET_MAX) stop(sprintf(
    'median log2FC = %+.3f: the contrast is not centred. Re-normalize (proteomics/quantification) before reading this table.',
    offset))
  offset
}

if (sys.nframe() == 0) {  # run as a script, not source()d
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 5) stop('Usage: Rscript centring_checks.R peptide_wide.csv samples.csv TRUE|FALSE Control Treatment')
  suppressPackageStartupMessages({ library(QFeatures); library(msqrob2) })
  peptide_wide <- read.csv(args[1], check.names = FALSE, stringsAsFactors = FALSE)
  sample_info <- read.csv(args[2], stringsAsFactors = FALSE)
  do_norm <- as.logical(args[3])
  ref_level <- args[4]; test_level <- args[5]
  runs <- sample_info$run
  col_data <- data.frame(quantCols = runs, condition = factor(sample_info$condition, levels = c(ref_level, test_level)),
                         sample = factor(runs), row.names = runs)
  pe <- readQFeatures(assayData = peptide_wide, quantCols = runs, colData = col_data, name = 'peptideRaw', verbose = FALSE)
  pe <- zeroIsNA(pe, 'peptideRaw')
  pe <- logTransform(pe, base = 2, i = 'peptideRaw', name = 'peptideLog')
  rowData(pe[['peptideLog']])$nNonZero <- rowSums(!is.na(assay(pe[['peptideLog']])))
  pe <- filterFeatures(pe, ~ nNonZero >= 2, keep = TRUE)
  cond <- colData(pe)$condition
  src <- 'peptideLog'
  if (do_norm) {
    pe <- normalize(pe, i = 'peptideLog', name = 'peptideNorm', method = 'center.median')
    sd_ratio <- check_sd_ratio(assay(pe[['peptideLog']]), assay(pe[['peptideNorm']]), cond)
    cat(sprintf('residual SD ratio (after / before) = %.3f\n', sd_ratio))
    src <- 'peptideNorm'
  }
  pe <- suppressWarnings(aggregateFeatures(pe, i = src, fcol = 'protein', name = 'protein',
                                           fun = MsCoreUtils::robustSummary, na.rm = TRUE))
  pe <- suppressWarnings(msqrob(pe, i = 'protein', formula = ~condition, robust = TRUE))
  contrast_name <- paste0('condition', test_level)
  L <- makeContrast(paste(contrast_name, '= 0'), parameterNames = contrast_name)
  pe <- hypothesisTest(pe, i = 'protein', contrast = L)
  res <- rowData(pe[['protein']])[[contrast_name]]
  res <- res[!is.na(res$adjPval), ]
  cat(sprintf('tested %d | calls (adjPval < 0.05) %d | median logFC %+.4f\n',
              nrow(res), sum(res$adjPval < 0.05), median(res$logFC)))
  offset <- check_offset(res$logFC)
  cat(sprintf('offset check passed (|%+.4f| <= %.2f)\n', offset, OFFSET_MAX))
}
