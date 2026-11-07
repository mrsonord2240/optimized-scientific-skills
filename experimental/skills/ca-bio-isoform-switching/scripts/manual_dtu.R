args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) stop('Usage: Rscript manual_dtu.R salmon_dir metadata.tsv tx2gene.tsv results.csv')
salmon_dir <- args[1]
meta <- read.delim(args[2], stringsAsFactors=FALSE, check.names=FALSE)
tx2gene <- read.delim(args[3], stringsAsFactors=FALSE, check.names=FALSE)
if (!all(c('sample_id','condition') %in% names(meta)) || anyDuplicated(meta$sample_id) || anyNA(meta$sample_id) || anyNA(meta$condition)) stop('metadata needs unique sample_id and condition')
if (length(unique(meta$condition)) != 2 || min(table(meta$condition)) < 2) stop('need two conditions and at least two replicates each')
if (!all(c('tx','gene') %in% names(tx2gene)) || anyDuplicated(tx2gene$tx)) stop('tx2gene needs unique tx and gene columns')
files <- setNames(file.path(salmon_dir, meta$sample_id, 'quant.sf'), meta$sample_id)
if (any(!file.exists(files))) stop('missing quant.sf input')
library(tximport); library(DRIMSeq); library(DEXSeq); library(stageR)

# meta: data.frame(sample_id, condition); files: named vector of quant.sf paths (names = sample_id)
# tx2gene: data.frame(tx = transcript ID, gene = gene ID)
# DRIMSeq normalizes count-column names. Normalize both inputs here so IDs such as
# `1-ctrl` remain aligned; keep the original metadata separately for display/reporting.
stopifnot(identical(names(files), meta$sample_id))
manual_sample_id <- make.names(meta$sample_id, unique = TRUE)
names(files) <- meta$sample_id <- manual_sample_id
txi <- tximport(files, type = 'salmon', txOut = TRUE, countsFromAbundance = 'no')
cts <- txi$counts
txdf <- data.frame(gene_id = tx2gene$gene[match(rownames(cts), tx2gene$tx)], feature_id = rownames(cts), cts, check.names = FALSE)
txdf <- txdf[!is.na(txdf$gene_id), ]

samples <- data.frame(sample_id = colnames(cts), condition = factor(meta$condition[match(colnames(cts), meta$sample_id)]))
stopifnot(!anyNA(samples$condition))

n <- nrow(samples); n_small <- min(table(samples$condition))
d <- dmDSdata(counts = txdf, samples = samples)
d <- dmFilter(d, min_samps_gene_expr = n, min_gene_expr = 10,
              min_samps_feature_expr = n_small, min_feature_expr = 10,
              min_samps_feature_prop = n_small, min_feature_prop = 0.1)

# Qualify with DRIMSeq:: -- once DEXSeq is attached, Biobase::samples masks DRIMSeq::samples() ("unable to find an
# inherited method for function 'samples' for signature 'object = "dmDSdata"'"); counts() is also exported by
# DEXSeq, DESeq2 and BiocGenerics, so qualify it too.
cnt <- DRIMSeq::counts(d)
dxd <- DEXSeqDataSet(
    countData = round(as.matrix(cnt[, samples$sample_id])),
    sampleData = DRIMSeq::samples(d),
    design = ~ sample + exon + condition:exon,
    featureID = cnt$feature_id,
    groupID = cnt$gene_id
)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, quiet = TRUE)
dxd <- testForDEU(dxd, reducedModel = ~ sample + exon)
dxr <- DEXSeqResults(dxd, independentFiltering = FALSE)
qval <- perGeneQValue(dxr)

pConfirmation <- matrix(dxr$pvalue, ncol = 1, dimnames = list(dxr$featureID, NULL))
pConfirmation[is.na(pConfirmation)] <- 1
tx2gene_d <- data.frame(transcript = dxr$featureID, gene = dxr$groupID)

stageRObj <- stageRTx(pScreen = qval, pConfirmation = pConfirmation, pScreenAdjusted = TRUE, tx2gene = tx2gene_d)
stageRObj <- stageWiseAdjustment(stageRObj, method = 'dtu', alpha = 0.05)

results <- getAdjustedPValues(stageRObj, order = FALSE, onlySignificantGenes = FALSE)   # geneID, txID, gene, transcript

write.csv(results, args[4], row.names=FALSE)
