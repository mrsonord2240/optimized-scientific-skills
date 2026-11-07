#!/usr/bin/env Rscript
# Purpose: msqrob2 differential abundance from a MaxQuant evidence.txt: drop flagged rows, build a QFeatures
#   object, log2, robust summarization to protein (or a peptide-level mixed model), msqrob, contrast test.
# Inputs:  evidence.txt; annotation CSV (Raw.file, Condition, BioReplicate, IsotopeLabelType); contrast
#          Treatment-Control (condition names from the annotation); output CSV.
# Usage:   Rscript msqrob2_de.R evidence.txt annotation.csv Treatment-Control results.csv [model=summary|peptide] [offset_max=0.05]
#            model=summary  robustSummary to protein, then msqrob (robust = TRUE); default
#            model=peptide  every peptide kept as a degree of freedom: ~condition + (1|sample) + (1|feature)
# Output:  results CSV (protein, logFC, se, df, t, pval, adjPval) and <results>_undetected.csv (proteins with no
#          observation in a condition, and proteins msqrob2 could not fit: adjPval NA). Stops if the median
#          logFC is above offset_max. The script does not normalize the peptide table.
# Checked: not executed. Written against msqrob2 1.14.1, QFeatures 1.16.0, MsCoreUtils 1.18.0 (R 4.4.3 / Bioconductor 3.20).
suppressPackageStartupMessages({ library(QFeatures); library(msqrob2) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop('Usage: Rscript msqrob2_de.R evidence.txt annotation.csv Treatment-Control results.csv [model=summary|peptide] [offset_max=0.05]')
opt <- list(model = 'summary', offset_max = '0.05')
for (a in args[-(1:4)]) { i <- regexpr('=', a, fixed = TRUE); opt[[substr(a, 1, i - 1)]] <- substr(a, i + 1, nchar(a)) }
parts <- strsplit(args[3], '-', fixed = TRUE)[[1]]
test_level <- parts[1]; ref_level <- parts[2]

ann <- read.csv(args[2])
if (!all(c(test_level, ref_level) %in% ann$Condition)) stop('Contrast names must be values of the Condition column')
runs <- as.character(ann$Raw.file)
ev <- read.table(args[1], sep = '\t', header = TRUE, quote = '', comment.char = '')
# Empty flag cells read as logical NA and NA != '+' is NA, which would drop every row: use %in%.
ev <- ev[!(ev$Reverse %in% '+') & !(ev$Potential.contaminant %in% '+') & !is.na(ev$Intensity) & ev$Intensity > 0, ]
ev$feature <- paste(ev$Modified.sequence, ev$Charge, sep = '_')
agg <- aggregate(Intensity ~ feature + Raw.file + Leading.razor.protein, data = ev, FUN = sum)
wide <- reshape(agg, idvar = c('feature', 'Leading.razor.protein'), timevar = 'Raw.file', direction = 'wide')
colnames(wide) <- sub('Intensity.', '', colnames(wide), fixed = TRUE)
wide <- wide[, c('feature', 'Leading.razor.protein', runs)]
names(wide)[2] <- 'protein'
cat('peptides:', nrow(wide), '| runs:', length(runs), '| proteins:', length(unique(wide$protein)), '\n')

col_data <- data.frame(quantCols = runs, condition = factor(ann$Condition, levels = c(ref_level, test_level)),
                       sample = factor(runs), row.names = runs)
pe <- readQFeatures(assayData = wide, quantCols = runs, colData = col_data, name = 'peptideRaw', verbose = FALSE)
pe <- zeroIsNA(pe, 'peptideRaw')
pe <- logTransform(pe, base = 2, i = 'peptideRaw', name = 'peptideLog')
rowData(pe[['peptideLog']])$nNonZero <- rowSums(!is.na(assay(pe[['peptideLog']])))
pe <- filterFeatures(pe, ~ nNonZero >= 2, keep = TRUE)

cond <- colData(pe)$condition  # colData lives on the QFeatures object, not on the assay
obs <- sapply(levels(cond), function(g)
  tapply(rowSums(!is.na(assay(pe[['peptideLog']])[, cond == g, drop = FALSE])), rowData(pe[['peptideLog']])$protein, sum))
undetected <- rownames(obs)[apply(obs, 1, min) == 0]

contrast_name <- paste0('condition', test_level)
L <- makeContrast(paste(contrast_name, '= 0'), parameterNames = contrast_name)
if (opt$model == 'peptide') {
  pe <- suppressWarnings(msqrobAggregate(pe, i = 'peptideLog', fcol = 'protein', name = 'proteinLmer',
                         formula = ~condition + (1 | sample) + (1 | feature), ridge = FALSE))  # ridge = TRUE is refused on a two-group design
  assay_name <- 'proteinLmer'
} else {
  pe <- suppressWarnings(aggregateFeatures(pe, i = 'peptideLog', fcol = 'protein', name = 'protein',
                                           fun = MsCoreUtils::robustSummary, na.rm = TRUE))
  pe <- suppressWarnings(msqrob(pe, i = 'protein', formula = ~condition, robust = TRUE))
  assay_name <- 'protein'
}
pe <- hypothesisTest(pe, i = assay_name, contrast = L)
res <- rowData(pe[[assay_name]])[[contrast_name]]
res <- data.frame(protein = rownames(pe[[assay_name]]), as.data.frame(res), row.names = NULL)
untestable <- res$protein[is.na(res$adjPval)]
res <- res[!is.na(res$adjPval), ]
res <- res[order(res$pval), ]

offset <- median(res$logFC, na.rm = TRUE)
cat(sprintf('tested %d | undetected in one condition %d | median logFC %+.4f | significant (adjPval < 0.05) %d\n',
            nrow(res), length(undetected), offset, sum(res$adjPval < 0.05)))
if (abs(offset) > as.numeric(opt$offset_max))
  stop(sprintf('median logFC = %+.3f: the contrast is not centred. Re-normalize the peptide table on features present in every run, or at protein level, before reading this table.', offset))
write.csv(res, args[4], row.names = FALSE)
write.csv(data.frame(protein = c(undetected, setdiff(untestable, undetected)),
                     reason = c(rep('no observation in one condition', length(undetected)),
                                rep('not fittable (adjPval NA)', length(setdiff(untestable, undetected))))),
          sub('\\.csv$', '_undetected.csv', args[4]), row.names = FALSE)
