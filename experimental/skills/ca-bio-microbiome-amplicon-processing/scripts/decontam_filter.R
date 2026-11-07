#!/usr/bin/env Rscript
# Checked on decontam (Bioconductor 3.18+)
# Remove reagent/kit contaminant ASVs with decontam::isContaminant.
#
# Usage: Rscript decontam_filter.R seqtab_nochim.rds metadata.csv out_dir [threshold]
#   seqtab_nochim.rds  samples (rows) x ASVs (columns), as written by dada2_asv.R
#   metadata.csv       first column = sample name (matches rownames of the table);
#                      columns: is_control (TRUE/FALSE, negative controls) and/or dna_conc (qPCR/Qubit per sample)
#   threshold          default 0.1; 0.5 is aggressive
# Method follows the columns present: both -> combined; only is_control -> prevalence; only dna_conc -> frequency.
# Writes contaminants.csv and seqtab_clean.rds in out_dir.

suppressPackageStartupMessages(library(decontam))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript decontam_filter.R seqtab_nochim.rds metadata.csv out_dir [threshold]")
threshold <- if (length(args) >= 4) as.numeric(args[4]) else 0.1
dir.create(args[3], showWarnings = FALSE, recursive = TRUE)

seqtab <- readRDS(args[1])
meta <- read.csv(args[2], row.names = 1, stringsAsFactors = FALSE)
missing <- setdiff(rownames(seqtab), rownames(meta))
if (length(missing) > 0) stop("samples without metadata: ", paste(head(missing, 10), collapse = ", "))
meta <- meta[rownames(seqtab), , drop = FALSE]

has_neg <- "is_control" %in% names(meta)
has_conc <- "dna_conc" %in% names(meta)
if (!has_neg && !has_conc) stop("metadata needs is_control and/or dna_conc")

contam <- if (has_neg && has_conc) {
  isContaminant(seqtab, neg = as.logical(meta$is_control), conc = meta$dna_conc, method = "combined", threshold = threshold)
} else if (has_neg) {
  isContaminant(seqtab, neg = as.logical(meta$is_control), method = "prevalence", threshold = threshold)
} else {
  isContaminant(seqtab, conc = meta$dna_conc, method = "frequency", threshold = threshold)
}
contam$asv <- colnames(seqtab)
write.csv(contam, file.path(args[3], "contaminants.csv"), row.names = FALSE)

seqtab_clean <- seqtab[, !contam$contaminant, drop = FALSE]
saveRDS(seqtab_clean, file.path(args[3], "seqtab_clean.rds"))
cat(sprintf("ASVs removed: %d of %d; reads removed: %d of %d (%.2f%%)\n",
            sum(contam$contaminant), ncol(seqtab), sum(seqtab[, contam$contaminant]), sum(seqtab),
            100 * sum(seqtab[, contam$contaminant]) / sum(seqtab)))
cat("Done:", file.path(args[3], "seqtab_clean.rds"), "\n")
