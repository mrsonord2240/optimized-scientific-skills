# msstatstmt_multiplex.R -- multi-plex TMT from MaxQuant, bridged by each plex's reference channel.
# Purpose : the multi-plex route of bio-workflows-proteomics-pipeline (MSstatsTMT; never concatenate plexes).
# Inputs  : evidence.txt, proteinGroups.txt, annotation.csv (Run, Fraction, TechRepMixture, Channel,
#           Condition, Mixture, BioReplicate; the pooled reference channel is Condition = 'Norm')
# Usage   : Rscript msstatstmt_multiplex.R [evidence.txt] [proteinGroups.txt] [annotation.csv] [out.csv]
# Output  : out.csv = MSstatsTMT ComparisonResult plus adj.pvalue.global (BH across all contrasts)
args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) if (length(args) >= i) args[i] else default
evidence_file   <- arg(1, 'evidence.txt')
pg_file         <- arg(2, 'proteinGroups.txt')
annotation_file <- arg(3, 'annotation.csv')
out_file        <- arg(4, 'msstatstmt_results.csv')

library(MSstatsTMT)

# Multi-plex TMT from MaxQuant. Each plex (Mixture) carries a pooled reference channel, annotated
# Condition = 'Norm' (MSstatsTMT requires that exact label); that channel is the bridge. annotation.csv has one row per (Run, Channel) with
# columns Run, Fraction, TechRepMixture, Channel, Condition, Mixture, BioReplicate.
# quote = '' / comment.char = '' as in scripts/msstats_maxquant.R.
evidence <- read.table(evidence_file, sep = '\t', header = TRUE, quote = '', comment.char = '')
proteinGroups <- read.table(pg_file, sep = '\t', header = TRUE, quote = '', comment.char = '')
stopifnot(nrow(evidence) == length(readLines(evidence_file)) - 1,
          nrow(proteinGroups) == length(readLines(pg_file)) - 1)
annotation <- read.csv(annotation_file)
stopifnot('Norm' %in% annotation$Condition)   # no reference channel = nothing to bridge plexes with

tmt_input <- MaxQtoMSstatsTMTFormat(evidence, proteinGroups, annotation, use_log_file = FALSE,
                                     verbose = FALSE)

# Summarize to protein level. Within-plex global median normalization is followed by reference-channel
# normalization: every plex is rescaled to its own 'Norm' channel, which is the cross-plex (IRS-style)
# bridge, and the Norm channel is then dropped. Do not concatenate plexes before this step.
summ <- proteinSummarization(tmt_input, method = 'msstats', global_norm = TRUE, reference_norm = TRUE,
                             remove_norm_channel = TRUE, use_log_file = FALSE, verbose = FALSE)

# Contrasts: one row per comparison, columns = the Condition levels that survive (Norm is removed),
# in sorted order -- built from the levels, against the first one, as in scripts/limma_pipeline.R.
lv <- sort(setdiff(unique(as.character(annotation$Condition)), 'Norm'))
comparison <- t(sapply(lv[-1], function(l) as.numeric(lv == l) - as.numeric(lv == lv[1])))
colnames(comparison) <- lv
rownames(comparison) <- paste0(lv[-1], '_vs_', lv[1])
# moderated = TRUE borrows variance across proteins (limma-style); groupComparisonTMT adjusts within
# each contrast, so adjust ACROSS the rows yourself when you report several (as in scripts/msstats_maxquant.R).
tmt_res <- groupComparisonTMT(summ, contrast.matrix = comparison, moderated = TRUE,
                              adj.method = 'BH', use_log_file = FALSE, verbose = FALSE)$ComparisonResult
tmt_res$adj.pvalue.global <- p.adjust(tmt_res$pvalue, method = 'BH')
print(table(tmt_res$Label, tmt_res$adj.pvalue < 0.05))
write.csv(tmt_res, out_file, row.names = FALSE)
cat('wrote', out_file, '\n')
