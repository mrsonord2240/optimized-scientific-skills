args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop('usage: kegg_ora.R hits.txt universe.txt organism keyType result.csv')
library(clusterProfiler)
read_ids <- function(path) unique(trimws(readLines(path, warn = FALSE)))
hits <- read_ids(args[1]); hits <- hits[nzchar(hits)]
universe <- read_ids(args[2]); universe <- universe[nzchar(universe)]
if (!length(hits) || !length(universe) || !all(hits %in% universe)) stop('hits and universe must be nonempty and hits a subset of universe')
result <- enrichKEGG(gene = hits, organism = args[3], keyType = args[4], universe = universe,
                     pvalueCutoff = 0.05, pAdjustMethod = 'BH', minGSSize = 10,
                     maxGSSize = 500, qvalueCutoff = 0.2)
write.csv(as.data.frame(result), args[5], row.names = FALSE)
message('KEGG access date: ', Sys.Date(), '; input hits: ', length(hits), '; universe: ', length(universe))
