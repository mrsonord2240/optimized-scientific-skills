args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop('usage: consensus_csv.R counts.csv metadata.csv output.csv')
library(ALDEx2)
library(MicrobiomeStat)
counts <- as.matrix(read.csv(args[1], row.names = 1, check.names = FALSE))
meta <- read.csv(args[2], stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c('sample', 'Group') %in% names(meta)) || anyNA(meta[, c('sample', 'Group')]) ||
    anyDuplicated(meta$sample) || !setequal(colnames(counts), meta$sample)) stop('invalid sample metadata')
if (!is.numeric(counts) || anyNA(counts) || any(!is.finite(counts)) || any(counts < 0) ||
    any(counts != round(counts))) stop('counts must be nonnegative finite integers')
meta <- meta[match(colnames(counts), meta$sample), , drop = FALSE]
rownames(meta) <- meta$sample
if (length(unique(meta$Group)) != 2) stop('exactly two groups required')
meta$Group <- factor(meta$Group)
counts <- counts[rowSums(counts > 0) >= 0.10 * ncol(counts), , drop = FALSE]
if (!nrow(counts)) stop('no taxa pass prevalence filter')
set.seed(1)
aldex_out <- aldex(counts, as.character(meta$Group), mc.samples = 128, test = 't', effect = TRUE, denom = 'all')
linda_out <- linda(feature.dat = as.data.frame(counts), meta.dat = meta,
                   formula = '~ Group', feature.dat.type = 'count', prev.filter = 0, alpha = 0.05)
aldex_hits <- rownames(aldex_out)[!is.na(aldex_out$we.eBH) & !is.na(aldex_out$effect) & aldex_out$we.eBH < 0.05 & abs(aldex_out$effect) > 1]
linda_res <- linda_out$output[[1]]
linda_hits <- rownames(linda_res)[!is.na(linda_res$reject) & linda_res$reject]
taxa <- sort(union(aldex_hits, linda_hits))
write.csv(data.frame(taxon = taxa, in_aldex = taxa %in% aldex_hits,
                     in_linda = taxa %in% linda_hits), args[3], row.names = FALSE)
