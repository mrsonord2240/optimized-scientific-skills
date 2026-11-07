#!/usr/bin/env Rscript
# Purpose: condition DE within one cell type from single-cell counts: sum raw counts to one pseudobulk per
#   donor x condition, drop thin pseudobulks, fit DESeq2 with the donor as a blocking factor when donors
#   span both conditions, and write one table with gene symbols attached.
# Inputs:  10x directory (matrix.mtx[.gz] genes x cells, features.tsv[.gz], barcodes.tsv[.gz]) holding RAW
#          integer counts for ONE cell type (subset first; run once per cell type);
#          cell table (TSV/CSV; first column = barcode) with a donor column and a condition column.
# Usage:   Rscript pseudobulk_de.R <10x_dir> <cells.tsv> <out.tsv> [key=value ...]
#            donor=donor  condition=condition  ref=control   column names and the reference level
#            min_cells=10       drop a pseudobulk with fewer cells (Crowell 2020; OSCA multi-sample)
#            min_count=10       keep genes with at least this total raw count over the retained cells
#            all_genes_bh=false true = BH over every tested gene (independentFiltering and Cook's off)
#          or source() this file and call run_pseudobulk_de() to keep `dds`.
# Output:  gene_id, symbol, baseMean, log2FoldChange (unshrunken MLE), lfcSE, stat, pvalue, padj, and
#          mean_norm_<level> = mean size-factor-normalized pseudobulk count per condition, for a
#          caller-defined fold change such as log2(mean_norm_treated / mean_norm_control).
# Checked: DESeq2 1.46.0, Matrix 1.7.6 (R 4.4.1).
suppressPackageStartupMessages({ library(Matrix); library(DESeq2) })

read_10x <- function(dir) {
  pick <- function(stem) {
    hit <- file.path(dir, c(paste0(stem, '.gz'), stem))
    hit <- hit[file.exists(hit)]
    if (!length(hit)) stop('missing ', stem, '[.gz] in ', dir)
    hit[1]
  }
  counts <- as(readMM(pick('matrix.mtx')), 'CsparseMatrix')
  features <- read.delim(pick('features.tsv'), header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.delim(pick('barcodes.tsv'), header = FALSE, stringsAsFactors = FALSE)[[1]]
  if (nrow(counts) != nrow(features) || ncol(counts) != length(barcodes))
    stop('matrix is ', nrow(counts), ' x ', ncol(counts), ' but features/barcodes have ',
         nrow(features), ' / ', length(barcodes), ' rows')
  dimnames(counts) <- list(features[[1]], barcodes)
  symbols <- if (ncol(features) >= 2) features[[2]] else features[[1]]
  list(counts = counts, symbols = setNames(symbols, features[[1]]))
}

run_pseudobulk_de <- function(counts, symbols, cells, donor = 'donor', condition = 'condition',
                              ref = NULL, min_cells = 10, min_count = 10, all_genes_bh = FALSE) {
  for (column in c(donor, condition))
    if (!column %in% names(cells)) stop('cell table has no column "', column, '"')
  cells <- cells[match(colnames(counts), cells[[1]]), , drop = FALSE]
  if (anyNA(cells[[1]])) stop(sum(is.na(cells[[1]])), ' barcodes in the matrix are absent from the cell table')
  if (any(counts@x != round(counts@x))) stop('counts are not integers: pseudobulk needs RAW counts, not normalized values')

  cond <- factor(cells[[condition]])
  if (nlevels(cond) != 2) stop('expected two conditions, found: ', paste(levels(cond), collapse = ', '))
  if (!is.null(ref)) cond <- relevel(cond, ref = ref)  # first level is the reference
  unit <- paste(cells[[donor]], cond, sep = '|')
  n_cells <- table(unit)

  # Thin pseudobulks are mostly zeros and add variance; in a paired design the donor goes with them.
  thin <- names(n_cells)[n_cells < min_cells]
  keep_unit <- setdiff(names(n_cells), thin)
  donor_of <- sub('\\|[^|]*$', '', keep_unit)
  paired <- any(duplicated(donor_of))
  dropped_donors <- character(0)
  if (paired) {
    complete <- names(which(table(donor_of) == nlevels(cond)))
    dropped_donors <- setdiff(unique(as.character(cells[[donor]])), complete)
    keep_unit <- keep_unit[donor_of %in% complete]
  }
  keep_cell <- unit %in% keep_unit
  counts <- counts[, keep_cell, drop = FALSE]
  unit <- factor(unit[keep_cell])

  pb <- as.matrix(counts %*% sparse.model.matrix(~0 + unit))  # sum of raw counts per pseudobulk
  colnames(pb) <- levels(unit)
  storage.mode(pb) <- 'integer'
  coldata <- data.frame(
    donor = factor(sub('\\|[^|]*$', '', levels(unit))),
    condition = factor(sub('^.*\\|', '', levels(unit)), levels = levels(cond)),
    n_cells = as.integer(table(unit)), row.names = levels(unit))
  if (min(table(coldata$condition)) < 2) stop('fewer than 2 pseudobulks in a condition after filtering: no replication')

  tested <- rowSums(pb) >= min_count
  design <- if (paired) ~ donor + condition else ~ condition  # blocking factor first, condition last
  dds <- DESeqDataSetFromMatrix(pb[tested, , drop = FALSE], coldata, design)
  dds <- DESeq(dds, quiet = TRUE)
  coef <- paste0('condition_', levels(cond)[2], '_vs_', levels(cond)[1])
  res <- if (all_genes_bh) results(dds, name = coef, independentFiltering = FALSE, cooksCutoff = FALSE)
         else results(dds, name = coef)

  norm <- counts(dds, normalized = TRUE)
  means <- sapply(levels(cond), function(g) rowMeans(norm[, coldata$condition == g, drop = FALSE]))
  colnames(means) <- paste0('mean_norm_', levels(cond))
  table_out <- data.frame(gene_id = rownames(res), symbol = unname(symbols[rownames(res)]),
                          as.data.frame(res), means, check.names = FALSE, row.names = NULL)
  list(results = table_out[order(table_out$pvalue), ], dds = dds, coldata = coldata, paired = paired,
       coef = coef, dropped_pseudobulks = thin, dropped_donors = dropped_donors, n_tested = sum(tested))
}

if (sys.nframe() == 0) {  # run as a script, not source()d
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) stop('Usage: Rscript pseudobulk_de.R <10x_dir> <cells.tsv> <out.tsv> [donor=donor condition=condition ref=control min_cells=10 min_count=10 all_genes_bh=false]')
  opt <- list(donor = 'donor', condition = 'condition', ref = NULL, min_cells = '10', min_count = '10', all_genes_bh = 'false')
  for (pair in strsplit(args[-(1:3)], '=', fixed = TRUE)) {
    if (length(pair) != 2 || !pair[1] %in% names(opt)) stop('unknown option: ', paste(pair, collapse = '='))
    opt[[pair[1]]] <- pair[2]
  }
  tenx <- read_10x(args[1])
  cells <- read.delim(args[2], sep = if (grepl('\\.csv$', args[2])) ',' else '\t', stringsAsFactors = FALSE)
  out <- run_pseudobulk_de(tenx$counts, tenx$symbols, cells, opt$donor, opt$condition, opt$ref,
                           as.integer(opt$min_cells), as.integer(opt$min_count), tolower(opt$all_genes_bh) == 'true')
  write.table(out$results, args[3], sep = '\t', quote = FALSE, row.names = FALSE)
  cat('design:', if (out$paired) '~ donor + condition (paired)' else '~ condition (donors nested in condition)',
      '| coefficient:', out$coef, '\n')
  cat('pseudobulks kept:', nrow(out$coldata), '| cells per pseudobulk:', paste(range(out$coldata$n_cells), collapse = '-'),
      '| dropped (<', opt$min_cells, 'cells):', length(out$dropped_pseudobulks), '\n')
  if (length(out$dropped_donors)) cat('donors dropped (not in both conditions after the cell filter):', paste(out$dropped_donors, collapse = ', '), '\n')
  cat('genes tested:', out$n_tested, '| padj < 0.05:', sum(out$results$padj < 0.05, na.rm = TRUE), '| wrote', args[3], '\n')
}
