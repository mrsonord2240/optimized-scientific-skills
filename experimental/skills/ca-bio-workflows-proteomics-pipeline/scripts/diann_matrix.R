# diann_matrix.R -- DIA-NN report.parquet -> log2 protein-group matrix (PG.MaxLFQ), ready for limma.
# Purpose : the DIA-NN route of bio-workflows-proteomics-pipeline: 1% q-value filters at precursor, protein
#           group and global protein group BEFORE pivoting, then 0 -> NA, log2.
# Inputs  : report.parquet (DIA-NN 1.9+/2.x)
# Usage   : Rscript diann_matrix.R [report.parquet] [log2_matrix.csv]
# Output  : log2_matrix.csv (rows = Protein.Group, columns = Run); continue with limma on it as in SKILL.md.
args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) if (length(args) >= i) args[i] else default
report_file <- arg(1, 'report.parquet')
out_file    <- arg(2, 'diann_log2_matrix.csv')

library(arrow)
library(dplyr)
library(tidyr)

diann <- read_parquet(report_file)

# Filter to 1% FDR at precursor AND protein-group level before pivoting.
# Use the GLOBAL protein-group q-value for the cross-run matrix (per-run min(Q.Value) is anti-conservative).
# When MBR is ON, MBR has its own FDR: add the Lib.* q-values (Lib.Q.Value, Lib.PG.Q.Value <= 0.01).
diann_filt <- diann %>%
    filter(Q.Value <= 0.01 & PG.Q.Value <= 0.01 & Global.PG.Q.Value <= 0.01)

# PG.MaxLFQ is ALREADY cross-run MaxLFQ-normalized at report generation. Re-normalizing it
# double-normalizes -- go straight to log2 + limma with no further normalization. To apply the
# skill's own median-centering instead, pivot raw PG.Quantity here, not PG.MaxLFQ.
protein_matrix <- diann_filt %>%
    select(Protein.Group, Run, PG.MaxLFQ) %>%
    distinct() %>%
    pivot_wider(names_from = Run, values_from = PG.MaxLFQ)

# PG.MaxLFQ path: log2-transform and go straight to limma (no re-normalization). DIA-NN writes 0
# for "not quantified in this run", so 0 -> NA FIRST: log2(0) is -Inf, and -Inf cells propagate
# silently until eBayes stops with "missing value where TRUE/FALSE needed".
m <- as.matrix(protein_matrix[, -1])
rownames(m) <- protein_matrix$Protein.Group
m[m == 0] <- NA
log2_matrix <- log2(m)
stopifnot(!any(is.infinite(log2_matrix)))
write.csv(log2_matrix, out_file)
cat('wrote', out_file, ':', nrow(log2_matrix), 'protein groups x', ncol(log2_matrix), 'runs\n')
