#!/usr/bin/env Rscript
# Purpose: MSstats feature-level two-group test from a MaxQuant evidence.txt: MaxQtoMSstatsFormat, dataProcess
#   (TMP summary, no AFT imputation), groupComparison, with undetected (oneConditionMissing) proteins split out.
# Inputs:  evidence.txt, proteinGroups.txt, annotation CSV (Raw.file, Condition, BioReplicate, IsotopeLabelType),
#          reference and test condition names, normalization ('equalizeMedians' = MSstats default, or FALSE),
#          output prefix. `normalization` is the per-run median normalization that references/feature_level.md
#          warns about: run its centring checks, or pass FALSE and normalize at the protein level.
# Usage:   Rscript msstats_group_comparison.R evidence.txt proteinGroups.txt annotation.csv Control Treatment \
#            equalizeMedians out_prefix
#          writes out_prefix_tested.csv and out_prefix_undetected.csv.
# Checked: MSstats 4.14.2 (R 4.4.3 / Bioconductor 3.20).
suppressPackageStartupMessages(library(MSstats))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) stop('Usage: Rscript msstats_group_comparison.R evidence.txt proteinGroups.txt annotation.csv Control Treatment equalizeMedians|FALSE out_prefix')
evidence <- read.table(args[1], sep = '\t', header = TRUE, quote = '', comment.char = '')
protein_groups <- read.table(args[2], sep = '\t', header = TRUE, quote = '', comment.char = '')
annotation <- read.csv(args[3])
ref_level <- args[4]; test_level <- args[5]
normalization <- if (args[6] == 'FALSE') FALSE else args[6]

input <- MaxQtoMSstatsFormat(evidence = evidence, proteinGroups = protein_groups,
                             annotation = annotation, use_log_file = FALSE)  # annotation: Raw.file, Condition, BioReplicate, IsotopeLabelType
proc <- dataProcess(input,
                    normalization = normalization,  # per-run median normalization: run the centring checks; FALSE = none (audit set: -0.184 / 20.8% FDR vs +0.008 / 0.0%)
                    summaryMethod = 'TMP', censoredInt = 'NA', MBimpute = FALSE, use_log_file = FALSE)

lv <- levels(proc$ProteinLevelData$GROUP)
stopifnot(all(c(ref_level, test_level) %in% lv))
contrast <- matrix(ifelse(lv == test_level, 1, ifelse(lv == ref_level, -1, 0)), nrow = 1,
                   dimnames = list(paste0(test_level, '-', ref_level), lv))  # order = levels(GROUP)
res <- groupComparison(contrast.matrix = contrast, data = proc, use_log_file = FALSE)$ComparisonResult
res$Protein <- as.character(res$Protein)

undetected <- res[res$issue %in% 'oneConditionMissing', ]   # infinite log2FC; report as undetected, not as a ratio
tested <- res[is.finite(res$log2FC) & !is.na(res$adj.pvalue), ]
# columns: Protein, Label, log2FC, SE, Tvalue, DF, pvalue, adj.pvalue, issue (adj.pvalue is the BH p)
write.csv(tested, paste0(args[7], '_tested.csv'), row.names = FALSE)
write.csv(undetected, paste0(args[7], '_undetected.csv'), row.names = FALSE)
cat(sprintf('tested %d | undetected %d | calls (adj.pvalue < 0.05) %d | median log2FC %+.3f | median SE %.3f\n',
            nrow(tested), nrow(undetected), sum(tested$adj.pvalue < 0.05),
            median(tested$log2FC), median(tested$SE)))
