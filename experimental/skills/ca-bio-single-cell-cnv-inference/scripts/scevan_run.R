args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || length(args) > 4) stop("Usage: Rscript scevan_run.R counts.tsv sample output.rds [normals.txt]")
library(SCEVAN)
counts <- as.matrix(read.table(args[1], header = TRUE, row.names = 1, sep = "\t", check.names = FALSE))
opts <- list(counts, sample = args[2], par_cores = 4, SUBCLONES = TRUE,
             ClonalCN = TRUE, plotTree = FALSE, organism = "human", ngenes_chr = 5)
if (length(args) == 4) opts$norm_cell <- scan(args[4], what = character(), quiet = TRUE)
result <- do.call(pipelineCNA, opts)
saveRDS(result, args[3])
