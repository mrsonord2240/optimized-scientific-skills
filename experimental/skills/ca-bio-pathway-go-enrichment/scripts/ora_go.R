#!/usr/bin/env Rscript
# Purpose: GO over-representation analysis (enrichGO, or enricher on a custom TERM2GENE) of a gene LIST
#   against an explicit background universe; one results table.
# Usage:   Rscript ora_go.R genes.txt universe.txt results.csv [key=value ...]
#          Rscript ora_go.R de=de_results.csv results.csv [key=value ...]
#            genes.txt / universe.txt   one gene ID per line (foreground, and the genes that could have been selected)
#            de=<csv>                   DE table; foreground = padj < padj & |log2FoldChange| > lfc,
#                                       universe = rows with a non-NA pvalue. Columns: gene_col, padj, log2FoldChange, pvalue
#            gene_col=gene_id  direction=both|up|down   (with de=)
#            padj=0.05  lfc=1                           (with de=)
#            keytype=SYMBOL    SYMBOL|ENTREZID|ENSEMBL or any key of the OrgDb
#            orgdb=org.Hs.eg.db   (org.Mm.eg.db mouse, org.Dr.eg.db zebrafish, org.Sc.sgd.db yeast with keytype=ORF)
#            ont=BP            BP|MF|CC|ALL   (ALL runs the three separately and simplifies each)
#            simplify=true     collapse redundant GO terms, cutoff 0.7, per ontology
#            term2gene=<tsv>   two columns term, gene (same ID type as keytype): runs enricher, ignores ont/simplify
#            minGSSize=10  maxGSSize=500  pvalueCutoff=0.05  qvalueCutoff=0.2
# Output:  results.csv: ID, Description, GeneRatio, BgRatio, FoldEnrichment, pvalue, p.adjust, qvalue, geneID, Count
#          [, ONTOLOGY]; results.csv.log holds the foreground/universe sizes and the ID-conversion rate.
# Checked: clusterProfiler 4.18.4, org.Hs.eg.db 3.22 (not yet executed after the recut).
suppressPackageStartupMessages({library(clusterProfiler)})

args <- commandArgs(trailingOnly = TRUE)
opts <- list(gene_col = 'gene_id', direction = 'both', padj = '0.05', lfc = '1', keytype = 'SYMBOL',
             orgdb = 'org.Hs.eg.db', ont = 'BP', simplify = 'true', term2gene = '', de = '',
             minGSSize = '10', maxGSSize = '500', pvalueCutoff = '0.05', qvalueCutoff = '0.2')
positional <- character(0)
for (a in args) {
  kv <- strsplit(a, '=', fixed = TRUE)[[1]]
  if (length(kv) >= 2 && kv[1] %in% names(opts)) opts[[kv[1]]] <- paste(kv[-1], collapse = '=')
  else if (grepl('=', a, fixed = TRUE)) stop('unknown argument ', a)
  else positional <- c(positional, a)
}
use_de <- nzchar(opts$de)
if (length(positional) != (if (use_de) 1 else 3)) stop('usage: ora_go.R genes.txt universe.txt results.csv [key=value ...]  or  ora_go.R de=de.csv results.csv [key=value ...]')
out <- positional[length(positional)]

orgdb_name <- opts$orgdb
if (!requireNamespace(orgdb_name, quietly = TRUE)) stop('Install ', orgdb_name, ' first')
OrgDb <- getExportedValue(orgdb_name, orgdb_name)
if (!opts$keytype %in% AnnotationDbi::keytypes(OrgDb))
  stop('keytype=', opts$keytype, ' is not a key of ', orgdb_name, ': ', paste(AnnotationDbi::keytypes(OrgDb), collapse = ', '))
ont <- toupper(opts$ont)
if (!ont %in% c('BP', 'MF', 'CC', 'ALL')) stop('ont must be BP, MF, CC or ALL')

if (use_de) {
  de <- read.csv(opts$de, check.names = FALSE)
  for (column in c(opts$gene_col, 'padj', 'log2FoldChange', 'pvalue'))
    if (!column %in% names(de)) stop('DE table has no column "', column, '"')
  hit <- !is.na(de$padj) & de$padj < as.numeric(opts$padj) & abs(de$log2FoldChange) > as.numeric(opts$lfc)
  if (opts$direction == 'up') hit <- hit & de$log2FoldChange > 0
  if (opts$direction == 'down') hit <- hit & de$log2FoldChange < 0
  fg_raw <- unique(as.character(de[[opts$gene_col]][hit]))
  bg_raw <- unique(as.character(de[[opts$gene_col]][!is.na(de$pvalue)]))
} else {
  fg_raw <- unique(trimws(readLines(positional[1])))
  bg_raw <- unique(trimws(readLines(positional[2])))
  fg_raw <- fg_raw[nzchar(fg_raw)]; bg_raw <- bg_raw[nzchar(bg_raw)]
}
if (length(fg_raw) == 0) stop('foreground is empty')
if (length(bg_raw) == 0) stop('universe is empty: ORA needs the genes that could have been selected')

# Map foreground and universe with the same call; enrichGO and the custom route both take ENTREZID.
to_entrez <- function(ids) {
  if (opts$keytype == 'ENTREZID' || nzchar(opts$term2gene)) return(ids)
  unique(suppressMessages(bitr(ids, fromType = opts$keytype, toType = 'ENTREZID', OrgDb = OrgDb))$ENTREZID)
}
fg <- to_entrez(fg_raw); bg <- to_entrez(bg_raw)
rate_fg <- length(fg) / length(fg_raw); rate_bg <- length(bg) / length(bg_raw)
note <- sprintf('foreground %d IDs -> %d mapped (%.1f%%); universe %d IDs -> %d mapped (%.1f%%)',
                length(fg_raw), length(fg), 100 * rate_fg, length(bg_raw), length(bg), 100 * rate_bg)
writeLines(note, paste0(out, '.log')); message(note)
if (rate_fg < 0.85 || rate_bg < 0.85) warning('More than 15% of genes failed ID conversion; check keytype and species.')
if (length(fg) == 0) stop('No foreground gene could be mapped: check keytype= and orgdb=')

size_args <- list(minGSSize = as.numeric(opts$minGSSize), maxGSSize = as.numeric(opts$maxGSSize),
                  pAdjustMethod = 'BH', pvalueCutoff = as.numeric(opts$pvalueCutoff),
                  qvalueCutoff = as.numeric(opts$qvalueCutoff))

if (nzchar(opts$term2gene)) {
  t2g <- read.delim(opts$term2gene, header = FALSE, stringsAsFactors = FALSE)[, 1:2]
  names(t2g) <- c('term', 'gene')
  res <- do.call(enricher, c(list(gene = fg, TERM2GENE = t2g, universe = bg), size_args))
  table_out <- as.data.frame(res)
} else {
  one <- function(o) {
    r <- do.call(enrichGO, c(list(gene = fg, universe = bg, OrgDb = OrgDb, keyType = 'ENTREZID', ont = o, readable = TRUE), size_args))
    if (!is.null(r) && nrow(as.data.frame(r)) > 0 && tolower(opts$simplify) == 'true')
      r <- simplify(r, cutoff = 0.7, by = 'p.adjust', select_fun = min, measure = 'Wang')
    d <- as.data.frame(r)
    if (nrow(d) > 0 && o != 'ALL' && !'ONTOLOGY' %in% names(d)) d$ONTOLOGY <- o
    d
  }
  table_out <- if (ont == 'ALL') do.call(rbind, lapply(c('BP', 'MF', 'CC'), one)) else one(ont)
}
if (nrow(table_out) == 0) message('No term passed the cutoffs. Rerun with pvalueCutoff=1 qvalueCutoff=1 to inspect; judge on p.adjust.')
write.csv(table_out, out, row.names = FALSE)
message('wrote ', out, ' (', nrow(table_out), ' terms)')
cat('orgdb', orgdb_name, as.character(packageVersion(orgdb_name)), 'GO.db', as.character(packageVersion('GO.db')), '\n',
    file = paste0(out, '.log'), append = TRUE)
