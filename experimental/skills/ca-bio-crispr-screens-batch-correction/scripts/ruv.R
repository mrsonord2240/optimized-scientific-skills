args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4 || length(args) > 5) stop("Usage: Rscript ruv.R counts.tsv ntcs.txt factors.tsv adjusted.tsv [k]")
library(RUVSeq)
counts_df <- read.delim(args[1], row.names = 1, check.names = FALSE)
ntc_sgrna_names <- readLines(args[2], warn = FALSE)
ntc_rownames <- rownames(counts_df)[rownames(counts_df) %in% ntc_sgrna_names]
stopifnot(length(ntc_rownames) > 0)
k <- if (length(args) == 5) as.integer(args[5]) else 2L
seqset <- newSeqExpressionSet(counts = as.matrix(counts_df))
ruv_corrected <- RUVg(seqset, cIdx = ntc_rownames, k = k)
write.table(pData(ruv_corrected), args[3], sep = "\t", quote = FALSE, col.names = NA)
write.table(normCounts(ruv_corrected), args[4], sep = "\t", quote = FALSE, col.names = NA)
