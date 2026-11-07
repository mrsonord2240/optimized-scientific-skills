#!/usr/bin/env Rscript
# Purpose: PTXQC report over a MaxQuant txt/ folder (per-metric scores in [0,1], heatmap PDF, mzQC, YAML).
# Usage:   Rscript scripts/ptxqc_report.R <maxquant_txt_dir> [files_written.txt]
# Output:  PTXQC writes its report files into the txt folder; the second argument lists what it added.
#          HTML only when Pandoc is reachable.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop('usage: ptxqc_report.R <maxquant_txt_dir> [files_written.txt]')
if (!dir.exists(args[1])) stop('not a directory: ', args[1])
if (!file.exists(file.path(args[1], 'summary.txt'))) stop('summary.txt missing: pass the MaxQuant txt/ folder')
suppressPackageStartupMessages(library(PTXQC))
if (!rmarkdown::pandoc_available()) message('Pandoc not found: HTML report skipped; PDF, mzQC and YAML are still written.')
before <- list.files(args[1])
createReport(txt_folder = args[1])
added <- setdiff(list.files(args[1]), before)
writeLines(added, if (length(args) >= 2) args[2] else file.path(args[1], 'ptxqc_files.txt'))
cat('PTXQC wrote:', paste(added, collapse = ', '), '\n')
