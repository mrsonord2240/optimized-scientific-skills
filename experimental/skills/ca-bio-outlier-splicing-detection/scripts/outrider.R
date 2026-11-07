library(OUTRIDER); library(BiocParallel)

args <- commandArgs(TRUE)
if (length(args) != 3) stop('Usage: outrider.R counts.tsv patient_id output.tsv')
countTable <- read.table(args[1], header=TRUE, row.names=1)
ods <- OutriderDataSet(countData = countTable)

ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)
# OUTRIDER 1.24.0: estimateBestQ() returns the number q; 1.28.1: returns the object, q in metadata
q_best <- estimateBestQ(ods)
if (is(q_best, 'OutriderDataSet')) { ods <- q_best; q_best <- metadata(ods)[['optimalEncDim']] }
q_best <- max(q_best, 2)    # OUTRIDER() requires q > 1; the optimal-hard-threshold estimate can be 1
bp <- if (.Platform$OS.type == 'windows') SerialParam() else MulticoreParam(8)
ods <- OUTRIDER(ods, q = q_best, BPPARAM = bp)

res <- results(ods, padjCutoff = 0.05, zScoreCutoff = 0)
patient_outliers <- res[res$sampleID == args[2], ]

write.table(patient_outliers, args[3], sep='\t', quote=FALSE, row.names=FALSE)
