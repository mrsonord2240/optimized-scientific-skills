args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: Rscript copykat_run.R counts.tsv sample output.rds")
library(copykat)
counts <- as.matrix(read.table(args[1], header = TRUE, row.names = 1, sep = "\t", check.names = FALSE))
res <- copykat(rawmat = counts, id.type = "S", ngene.chr = 5, win.size = 25,
  KS.cut = 0.1, sam.name = args[2], distance = "euclidean", norm.cell.names = "",
  genome = "hg20", n.cores = 4)
saveRDS(res, args[3])
