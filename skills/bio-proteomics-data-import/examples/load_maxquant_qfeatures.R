# Load a MaxQuant proteinGroups.txt into QFeatures at protein-group level.
# QFeatures makes row-data names syntactic, so filter with dotted MaxQuant flags.
suppressPackageStartupMessages(library(QFeatures))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop('Usage: Rscript load_maxquant_qfeatures.R proteinGroups.txt', call. = FALSE)
pg <- read.delim(args[[1]], quote = '', check.names = FALSE)
lfq <- grep('^LFQ intensity ', names(pg), value = TRUE)
if (!length(lfq)) stop('No LFQ intensity columns: this example is for label-free MaxQuant proteinGroups.txt', call. = FALSE)
required <- c('Reverse', 'Potential contaminant', 'Only identified by site')
if (!all(required %in% names(pg))) stop('Expected MaxQuant proteinGroups.txt bookkeeping columns', call. = FALSE)

qf <- readQFeatures(pg, quantCols = lfq, name = 'proteinGroups')
rd <- rowData(qf[[1]])
colnames(rd) <- make.names(colnames(rd))
rowData(qf[[1]]) <- rd
qf <- filterFeatures(qf, ~ !(Reverse %in% '+') & !(Potential.contaminant %in% '+') & !(Only.identified.by.site %in% '+'))
qf <- zeroIsNA(qf, 1)
qf <- logTransform(qf, i = 1, name = 'log2LFQ')
mat <- assay(qf[['log2LFQ']])
cat(sprintf('QFeatures log2 LFQ: %d protein groups x %d samples | -Inf: %s\n', nrow(mat), ncol(mat), any(is.infinite(mat))))
