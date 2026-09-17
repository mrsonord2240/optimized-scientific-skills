#!/usr/bin/env Rscript
# Feature-level differential abundance from a MaxQuant evidence.txt, with msqrob2.
# Checked on msqrob2 1.14.1, QFeatures 1.16.0, MsCoreUtils 1.18.0, limma 3.62.2 (R 4.4.3 / Bioc 3.20).
#
# Usage: Rscript msqrob2_peptide_level.R evidence.txt annotation.csv
#   annotation.csv: Raw.file, Condition, BioReplicate, IsotopeLabelType (the MSstats annotation layout)
# Falls back to a small simulated peptide table when no arguments are given.

suppressPackageStartupMessages({
  library(QFeatures)
  library(msqrob2)
})

CONTROL <- 'Control'
TREATMENT <- 'Treatment'

read_evidence <- function(evidence_path, annotation_path) {
  ev <- read.table(evidence_path, sep = '\t', header = TRUE, quote = '', comment.char = '')
  # Empty MaxQuant flag columns read as logical NA, and NA != '+' is NA, which would drop every row.
  ev <- ev[!(ev$Reverse %in% '+') & !(ev$Potential.contaminant %in% '+') &
             !is.na(ev$Intensity) & ev$Intensity > 0, ]
  ev$feature <- paste(ev$Modified.sequence, ev$Charge, sep = '_')
  ann <- read.csv(annotation_path)
  runs <- as.character(ann$Raw.file)
  agg <- aggregate(Intensity ~ feature + Raw.file + Leading.razor.protein, data = ev, FUN = sum)
  wide <- reshape(agg, idvar = c('feature', 'Leading.razor.protein'),
                  timevar = 'Raw.file', direction = 'wide')
  colnames(wide) <- sub('Intensity.', '', colnames(wide), fixed = TRUE)
  wide <- wide[, c('feature', 'Leading.razor.protein', runs)]
  names(wide)[2] <- 'protein'
  list(peptide_wide = wide, runs = runs,
       condition = factor(ann$Condition, levels = c(CONTROL, TREATMENT)))
}

simulate_peptides <- function(n_protein = 300, n_rep = 4, seed = 1) {
  set.seed(seed)
  runs <- c(paste0('C', seq_len(n_rep)), paste0('T', seq_len(n_rep)))
  proteins <- sprintf('P%04d', seq_len(n_protein))
  effect <- rep(0, n_protein)
  effect[1:20] <- 1.5
  effect[21:40] <- -1.5
  rows <- do.call(rbind, lapply(seq_len(n_protein), function(i) {
    n_pep <- sample(2:6, 1)
    base <- rnorm(1, 24, 2)
    do.call(rbind, lapply(seq_len(n_pep), function(j) {
      mu <- base + rnorm(1, 0, 1)
      y <- mu + c(rep(0, n_rep), rep(effect[i], n_rep)) + rnorm(2 * n_rep, 0, 0.3)
      data.frame(feature = sprintf('%s_pep%d', proteins[i], j), protein = proteins[i],
                 t(setNames(2^y, runs)), check.names = FALSE)
    }))
  }))
  list(peptide_wide = rows, runs = runs,
       condition = factor(rep(c(CONTROL, TREATMENT), each = n_rep), levels = c(CONTROL, TREATMENT)),
       truth = data.frame(protein = proteins, effect = effect))
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  d <- if (length(args) >= 2) read_evidence(args[1], args[2]) else simulate_peptides()
  runs <- d$runs
  cat('peptides:', nrow(d$peptide_wide), ' runs:', length(runs),
      ' proteins:', length(unique(d$peptide_wide$protein)), '\n')

  col_data <- data.frame(quantCols = runs, condition = d$condition,
                         sample = factor(runs), row.names = runs)
  pe <- readQFeatures(assayData = d$peptide_wide, quantCols = runs,
                      colData = col_data, name = 'peptideRaw', verbose = FALSE)
  pe <- zeroIsNA(pe, 'peptideRaw')
  pe <- logTransform(pe, base = 2, i = 'peptideRaw', name = 'peptideLog')
  rowData(pe[['peptideLog']])$nNonZero <- rowSums(!is.na(assay(pe[['peptideLog']])))
  pe <- filterFeatures(pe, ~ nNonZero >= 2, keep = TRUE)

  # Undetected proteins, taken before aggregation: msqrob2 returns adjPval = NA for them, not a list.
  cond <- colData(pe)$condition  # colData is on the QFeatures object, not on the individual assay
  obs <- sapply(levels(cond), function(g)
    tapply(rowSums(!is.na(assay(pe[['peptideLog']])[, cond == g, drop = FALSE])),
           rowData(pe[['peptideLog']])$protein, sum))
  undetected <- rownames(obs)[apply(obs, 1, min) == 0]
  cat('undetected in one condition:', length(undetected), '(reported separately, not as a fold change)\n')

  pe <- suppressWarnings(aggregateFeatures(pe, i = 'peptideLog', fcol = 'protein', name = 'protein',
                                           fun = MsCoreUtils::robustSummary, na.rm = TRUE))
  pe <- suppressWarnings(msqrob(pe, i = 'protein', formula = ~condition, robust = TRUE))
  contrast_name <- paste0('condition', TREATMENT)
  L <- makeContrast(paste(contrast_name, '= 0'), parameterNames = contrast_name)
  pe <- hypothesisTest(pe, i = 'protein', contrast = L)

  res <- rowData(pe[['protein']])[[contrast_name]]
  res$protein <- rownames(pe[['protein']])
  untestable <- res$protein[is.na(res$adjPval)]
  res <- res[!is.na(res$adjPval), ]
  cat('tested:', nrow(res), ' untestable (adjPval NA):', length(untestable), '\n')

  # A feature-level test has small enough SEs that a normalization offset becomes proteome-wide
  # significance. The median log2FC must be ~0 unless most of the proteome really moved.
  offset <- median(res$logFC, na.rm = TRUE)
  cat(sprintf('median logFC = %+.4f\n', offset))
  if (abs(offset) > 0.05)
    stop(sprintf('median logFC = %+.3f: the contrast is not centred. Re-normalize before reading this table.',
                 offset))

  sig <- res[res$adjPval < 0.05, ]
  sig <- sig[order(sig$pval), ]
  cat('significant at BH < 0.05:', nrow(sig), '\n')
  print(utils::head(as.data.frame(sig[, c('protein', 'logFC', 'se', 'df', 'adjPval')]), 5))

  if (!is.null(d$truth)) {
    fp <- sum(d$truth$effect[match(sig$protein, d$truth$protein)] == 0, na.rm = TRUE)
    cat(sprintf('planted truth: %d calls, %d false positives, realized FDR %.1f%%\n',
                nrow(sig), fp, ifelse(nrow(sig) > 0, 100 * fp / nrow(sig), 0)))
  }
  invisible(NULL)
}

main()
