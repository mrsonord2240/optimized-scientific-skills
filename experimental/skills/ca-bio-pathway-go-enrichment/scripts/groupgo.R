#!/usr/bin/env Rscript
# Purpose: descriptive GO-slim breakdown of a gene list at one DAG level (groupGO). Counts only: NOT a test.
# Usage:   Rscript groupgo.R genes.txt results.csv [keytype=ENTREZID] [orgdb=org.Hs.eg.db] [ont=BP] [level=3]
# Input:   genes.txt, one gene ID per line.
# Output:  results.csv: GO term, Description, Count, GeneRatio, geneID. No p-values exist.
# Checked: clusterProfiler 4.18.4 (not yet executed after the recut).
suppressPackageStartupMessages(library(clusterProfiler))

args <- commandArgs(trailingOnly = TRUE)
opts <- list(keytype = 'ENTREZID', orgdb = 'org.Hs.eg.db', ont = 'BP', level = '3')
positional <- character(0)
for (a in args) {
  kv <- strsplit(a, '=', fixed = TRUE)[[1]]
  if (length(kv) >= 2 && kv[1] %in% names(opts)) opts[[kv[1]]] <- paste(kv[-1], collapse = '=')
  else if (grepl('=', a, fixed = TRUE)) stop('unknown argument ', a)
  else positional <- c(positional, a)
}
if (length(positional) != 2) stop('usage: groupgo.R genes.txt results.csv [keytype=ENTREZID] [orgdb=org.Hs.eg.db] [ont=BP] [level=3]')
if (!requireNamespace(opts$orgdb, quietly = TRUE)) stop('Install ', opts$orgdb, ' first')
OrgDb <- getExportedValue(opts$orgdb, opts$orgdb)
genes <- unique(trimws(readLines(positional[1]))); genes <- genes[nzchar(genes)]
if (!length(genes)) stop('gene list is empty')
ggo <- groupGO(genes, OrgDb = OrgDb, keyType = opts$keytype, ont = toupper(opts$ont), level = as.integer(opts$level), readable = TRUE)
write.csv(as.data.frame(ggo), positional[2], row.names = FALSE)
message('wrote ', positional[2], ' (counts per GO term at level ', opts$level, '; no significance)')
