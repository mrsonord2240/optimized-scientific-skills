args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop('usage: reactome_ora.R hits.txt universe.txt organism output.csv')
library(ReactomePA)
read_ids <- function(path) unique(trimws(readLines(path, warn = FALSE)))
hits <- read_ids(args[1]); hits <- hits[nzchar(hits)]
universe <- read_ids(args[2]); universe <- universe[nzchar(universe)]
if (!length(hits) || !length(universe) || !all(hits %in% universe)) stop('hits and universe must be nonempty and hits a subset of universe')
result <- enrichPathway(gene = hits, organism = args[3], universe = universe,
                        pvalueCutoff = 0.05, qvalueCutoff = 0.2, pAdjustMethod = 'BH',
                        minGSSize = 10, maxGSSize = 500, readable = FALSE)
if (is.null(result)) stop('no mapped result; inspect identifiers and annotation')
write.csv(as.data.frame(result), args[4], row.names = FALSE)
message('Input hits: ', length(hits), '; universe: ', length(universe), '; result rows: ', nrow(as.data.frame(result)))
