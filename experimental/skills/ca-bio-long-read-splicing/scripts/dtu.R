args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4 || args[4] != 'gene_ids_without_underscores') stop('Usage: Rscript dtu.R reads_manifest.tsv flair_counts.tsv results.csv gene_ids_without_underscores')
library(DRIMSeq); library(stageR)   # do not attach DEXSeq here: it masks results()

manifest <- read.table(args[1], sep = '\t', col.names = c('sample_id', 'condition', 'batch', 'reads'), stringsAsFactors=FALSE)
counts <- read.table(args[2], header = TRUE, sep = '\t', check.names = FALSE)

if (anyDuplicated(manifest$sample_id) || anyNA(manifest$condition) || length(unique(manifest$condition)) != 2 || min(table(manifest$condition)) < 2) stop('manifest needs unique samples and two replicated conditions')

# First column = <isoform>_<gene> ('ids', or 'ID' with --sample_id_only); sample columns are <id>_<condition>_<batch> (just <id> with --sample_id_only)
sample_cols <- vapply(seq_len(nrow(manifest)), function(i) {
    hit <- intersect(c(manifest$sample_id[i], paste(manifest$sample_id[i], manifest$condition[i], manifest$batch[i], sep = '_')),
                     colnames(counts))
    stopifnot(length(hit) == 1)
    hit
}, character(1))
dm_counts <- data.frame(gene_id = sub('^.*_', '', counts[[1]]), feature_id = counts[[1]],
                        setNames(counts[, sample_cols], manifest$sample_id), check.names = FALSE)
samples <- data.frame(sample_id = manifest$sample_id, condition = factor(manifest$condition))

n <- nrow(samples); n_min <- min(table(samples$condition))
d <- dmDSdata(counts = dm_counts, samples = samples)
d <- dmFilter(d, min_samps_feature_expr = n_min, min_feature_expr = 5,
              min_samps_feature_prop = n_min, min_feature_prop = 0.1,
              min_samps_gene_expr = n, min_gene_expr = 10)

design <- model.matrix(~ condition, data = DRIMSeq::samples(d))
d <- dmPrecision(d, design = design)
d <- dmFit(d, design = design)
d <- dmTest(d, coef = colnames(design)[2])
res_gene <- DRIMSeq::results(d)
res_tx <- DRIMSeq::results(d, level = 'feature')

# stageR: gene-level screen, transcript-level confirmation
res_gene <- res_gene[!is.na(res_gene$pvalue), ]
res_tx <- res_tx[res_tx$gene_id %in% res_gene$gene_id & !is.na(res_tx$pvalue), ]
pScreen <- setNames(res_gene$pvalue, res_gene$gene_id)
pConfirm <- matrix(res_tx$pvalue, ncol = 1, dimnames = list(res_tx$feature_id, 'transcript'))
tx2gene <- data.frame(transcript = res_tx$feature_id, gene = res_tx$gene_id)
sr <- stageRTx(pScreen = pScreen, pConfirmation = pConfirm, pScreenAdjusted = FALSE, tx2gene = tx2gene)
sr <- stageWiseAdjustment(sr, method = 'dtu', alpha = 0.05, allowNA = TRUE)
dtu <- getAdjustedPValues(sr, order = TRUE, onlySignificantGenes = TRUE)   # stage-wise adjusted p per gene and transcript; NULL when no gene is significant

if (is.null(dtu)) dtu <- data.frame()
write.csv(dtu, args[3], row.names=FALSE)
