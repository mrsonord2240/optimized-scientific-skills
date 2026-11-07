# Summarize MaxQuant feature-level evidence to a normalized protein-level abundance matrix (MSstats TMP).
# Inputs : MaxQuant evidence.txt, proteinGroups.txt, and an MSstats annotation csv
#          (columns Raw.file, Condition, BioReplicate, IsotopeLabelType)
# Output : CSV of ProteinLevelData (one row per protein x run)
# Usage  : Rscript scripts/msstats_summarize.R evidence.txt proteinGroups.txt annotation.csv protein_abundance.csv [MBimpute=FALSE|TRUE]
# Checked: MSstats 4.14.2 (R 4.4.3)
suppressPackageStartupMessages(library(MSstats))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop('usage: msstats_summarize.R evidence.txt proteinGroups.txt annotation.csv out.csv [MBimpute=FALSE]')
mbimpute <- length(args) >= 5 && toupper(args[5]) == 'TRUE'

# quote = '' and comment.char = '': MaxQuant text fields contain apostrophes (5'-nucleotidase); default
# quoting silently truncates the table with only an 'EOF within quoted string' warning
evidence <- read.table(args[1], sep = '\t', header = TRUE, quote = '', comment.char = '')
protein_groups <- read.table(args[2], sep = '\t', header = TRUE, quote = '', comment.char = '')
stopifnot(nrow(evidence) == length(readLines(args[1])) - 1)  # every data line was read

maxquant_input <- MaxQtoMSstatsFormat(
    evidence = evidence,
    proteinGroups = protein_groups,
    annotation = read.csv(args[3])
)

# TMP = Tukey median polish. MBimpute = FALSE (default): no censoring model (censoredInt has no effect), and on/off
# proteins later come out of groupComparison as log2FC -Inf / issue 'oneConditionMissing'.
# MBimpute = TRUE is the AFT censored-imputation route (censoredInt = 'NA' marks NA intensities as censored).
processed <- dataProcess(maxquant_input, normalization = 'equalizeMedians',
                         summaryMethod = 'TMP', censoredInt = 'NA', MBimpute = mbimpute)

protein_abundance <- processed$ProteinLevelData
write.csv(protein_abundance, args[4], row.names = FALSE)
cat('proteins:', length(unique(protein_abundance$Protein)), '| rows:', nrow(protein_abundance), '\n')
