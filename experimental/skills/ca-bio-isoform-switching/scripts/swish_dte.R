args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop('Usage: Rscript swish_dte.R salmon_dir metadata.tsv results.csv')
salmon_dir <- args[1]
meta <- read.delim(args[2], stringsAsFactors=FALSE, check.names=FALSE)
if (!all(c('sample_id','condition') %in% names(meta)) || anyDuplicated(meta$sample_id) || anyNA(meta$sample_id) || anyNA(meta$condition)) stop('metadata needs unique sample_id and condition')
if (length(unique(meta$condition)) != 2 || min(table(meta$condition)) < 2) stop('need two conditions and at least two replicates each')
if (any(!file.exists(file.path(salmon_dir, meta$sample_id, 'quant.sf')))) stop('missing quant.sf input')
library(fishpond); library(tximeta)

# meta: data.frame(sample_id, condition); salmon_dir/<sample_id>/quant.sf from a run with --numGibbsSamples 20
coldata <- data.frame(names = meta$sample_id,
                      files = file.path(salmon_dir, meta$sample_id, 'quant.sf'),
                      condition = factor(meta$condition))   # swish needs a factor
se <- tximeta(coldata, skipMeta = TRUE)    # drop skipMeta once a linkedTxome matches your index
y <- scaleInfReps(se)
y <- labelKeep(y)
y <- y[mcols(y)$keep, ]

set.seed(1)
y <- swish(y, x = 'condition')
y <- computeInfRV(y)                        # adds mcols(y)$meanInfRV

dte_results <- as.data.frame(mcols(y))     # log2FC, pvalue, qvalue, meanInfRV
sig <- subset(dte_results, qvalue < 0.05)

dte_results$transcript_id <- rownames(dte_results)
write.csv(dte_results, args[3], row.names=FALSE)
