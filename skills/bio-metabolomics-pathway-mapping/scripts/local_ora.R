# Purpose: Local-Only ORA. Hypergeometric pathway enrichment computed locally from KEGGREST's public
#          pathway-to-compound table (generic reference data, not user data), against an explicit
#          assay-coverage background.
# Inputs:  hit KEGG IDs file (one per line, e.g. from map_compounds.R); reference (background) KEGG
#          IDs file (one per line); optional output CSV (default ora_local.csv); optional min hits (2).
# Output:  the CSV: pathway, total, hits, p.value, fdr (sorted by p.value).
# Usage:   Rscript scripts/local_ora.R kegg_ids.txt reference_metabolome.txt [ora_local.csv] [2]
# Note:    fetches keggLink('pathway','compound') from the live KEGG API (~3 s).
suppressMessages(library(KEGGREST))

local_kegg_ora <- function(hit_kegg_ids, universe_kegg_ids, min_hits = 2) {
  links <- keggLink('pathway', 'compound')             # public reference table; no user data sent
  cpd_ids  <- sub('^cpd:', '', names(links))
  path_ids <- sub('^path:map', 'hsa', unname(links))    # generic map#### -> organism-specific hsa####
  pw2cpd <- lapply(split(cpd_ids, path_ids), unique)

  universe_kegg_ids <- unique(universe_kegg_ids)
  pw2cpd <- lapply(pw2cpd, function(x) intersect(x, universe_kegg_ids))
  pw2cpd <- pw2cpd[lengths(pw2cpd) > 0]
  hits <- intersect(hit_kegg_ids, universe_kegg_ids)
  N <- length(universe_kegg_ids); k <- length(hits)

  out <- data.frame(
    pathway = names(pw2cpd),
    total   = lengths(pw2cpd),
    hits    = vapply(pw2cpd, function(s) length(intersect(s, hits)), integer(1))
  )
  out <- out[out$hits >= min_hits, ]
  out$p.value <- mapply(function(m, h) phyper(h - 1, m, N - m, k, lower.tail = FALSE), out$total, out$hits)
  out$fdr <- p.adjust(out$p.value, method = 'BH')
  out[order(out$p.value), ]
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop('usage: Rscript local_ora.R kegg_ids.txt reference_metabolome.txt [out.csv] [min_hits]')
kegg_ids <- readLines(args[1], warn = FALSE); kegg_ids <- kegg_ids[nzchar(kegg_ids)]
# reference_ids: the assay-coverage background -- one KEGG compound ID per line
reference_ids <- readLines(args[2], warn = FALSE); reference_ids <- reference_ids[nzchar(reference_ids)]
out_file <- if (length(args) >= 3) args[3] else 'ora_local.csv'
min_hits <- if (length(args) >= 4) as.integer(args[4]) else 2L

ora_local <- local_kegg_ora(kegg_ids, reference_ids, min_hits)
write.csv(ora_local, out_file, row.names = FALSE)
print(head(ora_local, 10))
