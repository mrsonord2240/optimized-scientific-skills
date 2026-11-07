#!/usr/bin/env Rscript
# Purpose: DESeq2 differential expression from a count matrix and a sample table, in one run: set the
#   reference level, build the design, fit, extract one named comparison, and write one results table.
# Inputs:  counts CSV (first column = gene ID, one column per sample, RAW integer counts);
#          samples CSV (first column = sample ID matching the count columns, then design columns).
# Usage:   Rscript deseq2_de.R counts.csv samples.csv results.csv ref=<level> [key=value ...]
#            condition=condition   column holding the groups to compare
#            ref=<level>           REQUIRED reference (denominator) level
#            versus=<level>        numerator level; required when the column has 3+ levels and test=wald
#            block=<col[,col]>     batch, subject or other covariates, placed before the condition
#            test=wald|lrt         lrt = "does this factor change the gene at all" across every level
#            shrink=none|apeglm|ashr   adds log2FoldChange_shrunk for ranking and plots
#            min_count=10  alpha=0.05  all_genes_bh=false
#          or source() this file and call run_deseq2() on your own DESeqDataSet (tximport input).
# Output:  gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj [, log2FoldChange_shrunk],
#          sorted by p-value. With test=lrt the fold-change columns are omitted: an LRT has none.
# Checked: DESeq2 1.46.0, apeglm 1.28.0 (R 4.4.1).
suppressPackageStartupMessages(library(DESeq2))

run_deseq2 <- function(dds, condition = 'condition', ref, versus = NULL, block = character(0),
                       test = c('wald', 'lrt'), shrink = c('none', 'apeglm', 'ashr'),
                       min_count = 10, alpha = 0.05, all_genes_bh = FALSE) {
  test <- match.arg(test); shrink <- match.arg(shrink)
  for (column in c(condition, block))
    if (!column %in% names(colData(dds))) stop('sample table has no column "', column, '"')
  group <- factor(colData(dds)[[condition]])
  if (!ref %in% levels(group)) stop('ref="', ref, '" is not a level of ', condition, ': ', paste(levels(group), collapse = ', '))
  colData(dds)[[condition]] <- relevel(group, ref = ref)  # set BEFORE DESeq(); the default is alphabetical
  for (column in block) if (!is.numeric(colData(dds)[[column]])) colData(dds)[[column]] <- factor(colData(dds)[[column]])
  design(dds) <- reformulate(c(block, condition))  # covariates first, the tested factor last
  dds <- dds[rowSums(counts(dds)) >= min_count, ]
  levels_now <- levels(colData(dds)[[condition]])

  if (test == 'lrt') {
    dds <- DESeq(dds, test = 'LRT', reduced = if (length(block)) reformulate(block) else ~ 1, quiet = TRUE)
    res <- results(dds, alpha = alpha)
    label <- paste('LRT: any difference among', paste(levels_now, collapse = ', '))
    table_out <- data.frame(gene = rownames(res), as.data.frame(res)[, c('baseMean', 'stat', 'pvalue', 'padj')])
  } else {
    if (is.null(versus)) {
      if (length(levels_now) != 2) stop(condition, ' has levels ', paste(levels_now, collapse = ', '), ': pass versus=<level>, or test=lrt')
      versus <- levels_now[2]
    }
    if (!versus %in% levels_now) stop('versus="', versus, '" is not a level of ', condition)
    dds <- DESeq(dds, quiet = TRUE)
    contrast <- c(condition, versus, ref)  # never results(dds) bare: it returns the LAST coefficient
    res <- if (all_genes_bh) results(dds, contrast = contrast, alpha = alpha, independentFiltering = FALSE, cooksCutoff = FALSE)
           else results(dds, contrast = contrast, alpha = alpha)
    label <- paste(versus, 'vs', ref)
    table_out <- data.frame(gene = rownames(res), as.data.frame(res))
    if (shrink != 'none') {  # shrunken LFC is for ranking and plots; pvalue and padj stay those of results()
      coef <- which(resultsNames(dds) == paste0(condition, '_', make.names(versus), '_vs_', make.names(ref)))
      shrunk <- if (shrink == 'apeglm' && length(coef) == 1) lfcShrink(dds, coef = coef, res = res, type = 'apeglm', quiet = TRUE)
                else lfcShrink(dds, contrast = contrast, res = res, type = 'ashr', quiet = TRUE)
      table_out$log2FoldChange_shrunk <- shrunk$log2FoldChange
    }
  }
  na_padj <- is.na(table_out$padj)
  list(results = table_out[order(table_out$pvalue), ], dds = dds, label = label,
       design = paste(deparse(design(dds)), collapse = ''), n_per_level = table(colData(dds)[[condition]]),
       n_sig = sum(table_out$padj < alpha, na.rm = TRUE),
       na_all_zero = sum(na_padj & table_out$baseMean == 0),
       na_outlier = sum(na_padj & table_out$baseMean > 0 & is.na(table_out$pvalue)),
       na_filtered = sum(na_padj & !is.na(table_out$pvalue)))
}

if (sys.nframe() == 0) {  # run as a script, not source()d
  args <- commandArgs(trailingOnly = TRUE)
  usage <- 'Usage: Rscript deseq2_de.R counts.csv samples.csv results.csv ref=<level> [condition=condition versus=<level> block=<col,col> test=wald|lrt shrink=none|apeglm|ashr min_count=10 alpha=0.05 all_genes_bh=false]'
  if (length(args) < 4) stop(usage)
  opt <- list(condition = 'condition', ref = NULL, versus = NULL, block = '', test = 'wald', shrink = 'none',
              min_count = '10', alpha = '0.05', all_genes_bh = 'false')
  for (pair in strsplit(args[-(1:3)], '=', fixed = TRUE)) {
    if (length(pair) != 2 || !pair[1] %in% names(opt)) stop('unknown option: ', paste(pair, collapse = '='), '\n', usage)
    opt[[pair[1]]] <- pair[2]
  }
  if (is.null(opt$ref)) stop('ref=<level> is required: DESeq2 would otherwise pick the alphabetically first level')
  counts <- as.matrix(read.csv(args[1], row.names = 1, check.names = FALSE))
  samples <- read.csv(args[2], row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  if (!setequal(colnames(counts), rownames(samples))) stop('sample IDs differ between the count columns and the sample table')
  if (any(counts != round(counts))) stop('counts are not integers: use raw counts, or tximport for Salmon/kallisto/RSEM (routes/tximport.md)')
  storage.mode(counts) <- 'integer'
  dds <- DESeqDataSetFromMatrix(counts[, rownames(samples), drop = FALSE], samples, ~ 1)
  out <- run_deseq2(dds, opt$condition, opt$ref, opt$versus, if (nzchar(opt$block)) strsplit(opt$block, ',')[[1]] else character(0),
                    opt$test, opt$shrink, as.integer(opt$min_count), as.numeric(opt$alpha), tolower(opt$all_genes_bh) == 'true')
  write.csv(out$results, args[3], row.names = FALSE)
  cat('design:', out$design, '| comparison:', out$label, '\n')
  cat('samples per level:', paste(names(out$n_per_level), out$n_per_level, sep = '=', collapse = ' '), '\n')
  cat('genes tested:', nrow(out$results), '| padj <', opt$alpha, ':', out$n_sig, '| padj NA: all-zero', out$na_all_zero,
      ', outlier', out$na_outlier, ', low-count filter', out$na_filtered, '\n')
  cat('wrote', args[3], '\n')
}
