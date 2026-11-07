#!/usr/bin/env Rscript
# Purpose: copy-number correction of a pooled CRISPR count table with CRISPRcleanR (unsupervised, no CN
#   profile needed), writing the corrected count table every hit caller must then read.
# Inputs:  raw guide count table (tab; column 1 = sgRNA, column 2 = gene, then samples; the control/baseline
#          sample is the FIRST sample column); a CRISPRcleanR library annotation (default: the built-in
#          KY_Library_v1.0, so pass library=<file> for any other library).
# Usage:   Rscript cn_correction.R <counts.txt> <out_prefix> [library=KY_Library_v1.0|<annotation.tsv>] [min_reads=30]
#            library file columns: CODE, GENES, EXONE, CHRM, STRAND, STARTpos, ENDpos (CRISPRcleanR format, tab-separated)
# Output:  <out_prefix>_cleanr_corrected_counts.txt  (same layout as the input, corrected counts)
# Checked: CRISPRcleanR 3.0.1, R 4.4.3 (A375 Project Score table, 3 replicates + plasmid; non-integer counts out).
suppressPackageStartupMessages(library(CRISPRcleanR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop('Usage: Rscript cn_correction.R <counts.txt> <out_prefix> [library=...] [min_reads=30]')
opt <- list(library = 'KY_Library_v1.0', min_reads = '30')
for (pair in strsplit(args[-(1:2)], '=', fixed = TRUE)) {
  if (length(pair) != 2 || !pair[1] %in% names(opt)) stop('unknown option: ', paste(pair, collapse = '='))
  opt[[pair[1]]] <- pair[2]
}
if (!file.exists(args[1])) stop('no such count file: ', args[1])
counts <- read.delim(args[1], check.names = FALSE, nrows = 2000)
if (ncol(counts) < 4) stop('need sgRNA, gene and at least two sample columns')
if (any(unlist(counts[-(1:2)]) != round(unlist(counts[-(1:2)])), na.rm = TRUE))
  stop('counts are not integers: give the raw .count.txt, not a normalized table')

if (opt$library == 'KY_Library_v1.0') {
  data(KY_Library_v1.0)
  library_annotation <- KY_Library_v1.0
} else {
  library_annotation <- read.delim(opt$library, row.names = 1, stringsAsFactors = FALSE)
}
if (!all(c('GENES', 'CHRM', 'STARTpos', 'ENDpos') %in% colnames(library_annotation)))
  stop('library annotation lacks CRISPRcleanR columns (GENES, CHRM, STARTpos, ENDpos)')

norm <- ccr.NormfoldChanges(args[1], min_reads = as.integer(opt$min_reads), EXPname = 'screen',
                            libraryAnnotation = library_annotation)   # arg 1 is the file PATH
gw_lfc <- ccr.logFCs2chromPos(norm$logFCs, library_annotation)        # $logFCs, not $norm_fold_changes
cleaned <- ccr.GWclean(gw_lfc, display = TRUE, label = args[2])
corrected <- ccr.correctCounts('screen', norm$norm_counts, cleaned, library_annotation)
# ccr.correctCounts returns an in-memory frame and writes nothing: persist it.
write.table(corrected, paste0(args[2], '_cleanr_corrected_counts.txt'), sep = '\t', quote = FALSE, row.names = FALSE)
cat('guides:', nrow(corrected), '| wrote', paste0(args[2], '_cleanr_corrected_counts.txt'), '\n')
