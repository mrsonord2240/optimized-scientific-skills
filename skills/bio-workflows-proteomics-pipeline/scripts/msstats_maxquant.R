# msstats_maxquant.R -- MaxQuant evidence.txt + proteinGroups.txt -> MSstats feature-level comparison.
# Purpose : the MSstats route of bio-workflows-proteomics-pipeline (two conditions, Control vs Treatment).
#           For three or more conditions replace the `comparison` matrix as the comment describes.
# Inputs  : evidence.txt, proteinGroups.txt, annotation.csv (Raw.file, Condition, BioReplicate, IsotopeLabelType)
# Usage   : Rscript msstats_maxquant.R [evidence.txt] [proteinGroups.txt] [annotation.csv] [out.csv]
# Output  : out.csv = MSstats ComparisonResult (Protein, Label, log2FC, pvalue, adj.pvalue, ...)
args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) if (length(args) >= i) args[i] else default
evidence_file   <- arg(1, 'evidence.txt')
pg_file         <- arg(2, 'proteinGroups.txt')
annotation_file <- arg(3, 'annotation.csv')
out_file        <- arg(4, 'msstats_comparison.csv')

library(MSstats)

# From MaxQuant. quote = '' and comment.char = '' are REQUIRED: MaxQuant text fields contain
# apostrophes (5'-nucleotidase) and '#', and the read.table defaults silently truncate the table
# with nothing but an 'EOF within quoted string' warning -- here 2954 of 10369 evidence rows and
# 62 of 296 proteins. Check the row count; do not trust the warning to stop you.
evidence <- read.table(evidence_file, sep = '\t', header = TRUE, quote = '', comment.char = '')
proteinGroups <- read.table(pg_file, sep = '\t', header = TRUE, quote = '', comment.char = '')
stopifnot(nrow(evidence) == length(readLines(evidence_file)) - 1,
          nrow(proteinGroups) == length(readLines(pg_file)) - 1)
annotation <- read.csv(annotation_file)

# Convert to MSstats format
msstats_input <- MaxQtoMSstatsFormat(evidence = evidence,
                                      proteinGroups = proteinGroups,
                                      annotation = annotation)

# Process data. 'equalizeMedians' carries the SAME symmetry assumption as the median centering in
# scripts/limma_pipeline.R: most proteins unchanged, and the changes roughly balanced up and down. On the
# 296-protein set above (44 up, 26 down) it left every true null shifted -0.19 log2 (t vs 0,
# p = 2e-55) and MSstats then correctly called 21 of the 179 true nulls at BH 5%, all negative.
# For a one-sided design use normalization = FALSE with an externally normalized input, or
# normalization = 'globalStandards' with globalStandardName = <spike-in / unchanged protein set>.
# Always check the null centre: mean log2FC over proteins you expect not to move should be ~0.
processed <- dataProcess(msstats_input, normalization = 'equalizeMedians',
                         summaryMethod = 'TMP', censoredInt = 'NA')

# Comparison. +1 on the numerator: Treatment=+1, Control=-1 so log2FC = Treatment - Control
# (positive = up in Treatment), matching the label and the limma makeContrasts(Treatment-Control) path.
# Columns must be ALL the Condition levels in sorted order, one ROW per contrast -- for three
# conditions sorted Control/HighDose/LowDose that is a 2 x 3 matrix, e.g.
#   rbind(HighDose_vs_Control = c(-1, 1, 0), LowDose_vs_Control = c(-1, 0, 1))
# with colnames c('Control','HighDose','LowDose'); groupComparison adjusts within each row, so
# adjust across the rows yourself (p.adjust on the pooled pvalue) when you report several.
comparison <- matrix(c(-1, 1), nrow = 1)
rownames(comparison) <- 'Treatment_vs_Control'
colnames(comparison) <- c('Control', 'Treatment')

results <- groupComparison(contrast.matrix = comparison, data = processed)
write.csv(results$ComparisonResult, out_file, row.names = FALSE)
cat('wrote', out_file, '\n')
