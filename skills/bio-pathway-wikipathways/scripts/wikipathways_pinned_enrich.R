# Purpose: reproducible WikiPathways ORA on a PINNED dated GMT (not the live current/ release).
#   Tries newest-first monthly releases inside the ~12-month retention window, splits the compound
#   name%version%wpid%org term, runs enricher(), prints the release date used (report it in methods).
# Inputs: sig_file = Entrez IDs of the significant genes, one per line; universe_file = Entrez IDs of
#   all tested genes, one per line; organism = exact WP name (default 'Homo sapiens'); out_csv optional.
# Usage: Rscript wikipathways_pinned_enrich.R sig_entrez.txt universe_entrez.txt ['Homo sapiens'] [out.csv]
suppressMessages({library(rWikiPathways); library(clusterProfiler); library(tidyr)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop('usage: Rscript wikipathways_pinned_enrich.R sig_entrez.txt universe_entrez.txt [organism] [out.csv]')
sig         <- unique(trimws(readLines(args[1])))
all_entrez  <- unique(trimws(readLines(args[2])))
organism    <- if (length(args) >= 3) args[3] else 'Homo sapiens'
out_csv     <- if (length(args) >= 4) args[4] else NULL
sig <- sig[nzchar(sig)]; all_entrez <- all_entrez[nzchar(all_entrez)]

# newest-first candidate release dates, all inside the ~12-month window (10th of each month)
candidates <- unique(format(Sys.Date() - seq(60, 330, by = 30), '%Y%m10'))
gmt <- NULL
for (archive_date in candidates) {   # a missing release 404s: step back one month, then fall back to Zenodo
  # downloadPathwayArchive needs an organism to actually download a file (organism=NULL opens the index)
  gmt <- tryCatch(suppressWarnings(downloadPathwayArchive(date = archive_date, organism = organism,
                                                          format = 'gmt', destpath = tempdir())),
                  error = function(e) NULL)
  if (!is.null(gmt) && file.exists(file.path(tempdir(), gmt))) break
  gmt <- NULL
}
if (is.null(gmt)) stop('no release in the last ~12 months; use the Zenodo GMT archive')

wp2gene <- read.gmt(file.path(tempdir(), gmt))
wp2gene <- separate(wp2gene, term, c('name', 'version', 'wpid', 'org'), sep = '%')   # term is a %-joined compound
t2g <- wp2gene[, c('wpid', 'gene')]   # TERM2GENE
t2n <- wp2gene[, c('wpid', 'name')]   # TERM2NAME

wp_pinned <- enricher(sig, universe = all_entrez, TERM2GENE = t2g, TERM2NAME = t2n)
cat('WikiPathways release used (report in methods):', archive_date, '\n')
res <- as.data.frame(wp_pinned)
print(head(res[, c('ID', 'Description', 'GeneRatio', 'BgRatio', 'p.adjust', 'Count')], 10))
if (!is.null(out_csv)) write.csv(res, out_csv, row.names = FALSE)
