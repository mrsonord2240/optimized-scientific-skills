# Purpose: size the number of DONORS for population-level DE from a scRNA-seq pilot:
#          sum each donor's cell counts per gene (pseudobulk), fit DESeq2 dispersions,
#          then ssizeRNA_vary on the per-gene mean/dispersion vectors.
# Inputs:  cell_counts RDS (named list, one genes x cells matrix per donor),
#          donor_condition CSV (columns donor, condition), then optional positional args:
#          fc (1.5), fdr (0.05), power (0.80), maxN (200), seed (20260918).
# Usage:   r.sh pseudobulk_donor_ssize.R cell_counts.rds donor_condition.csv [1.5] [0.05] [0.80] [200] [20260918]
# Checked: ssizeRNA 1.3.3, DESeq2 1.46.0, R 4.4.3 / Bioconductor 3.20.
suppressPackageStartupMessages({library(DESeq2); library(ssizeRNA)})
safe_ssize <- function(call) {
  tryCatch(call(), error = function(e) {
    if (identical(conditionMessage(e), "argument is of length zero")) return(NULL)
    stop(e)
  })
}
require_reachable_n <- function(res, maxN) {
  n <- if (is.null(res)) NA_real_ else res$ssize[, "ssize"]
  if (length(n) != 1L || is.na(n)) {
    stop(sprintf("no n <= %d reaches the target; raise maxN or revise fc/dispersion", maxN))
  }
  n
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: pseudobulk_donor_ssize.R cell_counts.rds donor_condition.csv [fc] [fdr] [power] [maxN] [seed]")
arg <- function(i, default) if (length(args) >= i) args[i] else default
fc    <- as.numeric(arg(3, 1.5)); fdr <- as.numeric(arg(4, 0.05))
power <- as.numeric(arg(5, 0.80)); maxN <- as.integer(arg(6, 200))
set.seed(as.integer(arg(7, 20260918)))

# cell_counts: named list, one genes x cells matrix per donor (from the scRNA-seq count matrix, split by donor)
cell_counts <- readRDS(args[1])
dc <- read.csv(args[2])
donor_condition <- dc$condition[match(names(cell_counts), dc$donor)]
pseudobulk <- sapply(cell_counts, rowSums)               # genes x donors -- sum, not mean, across cells
coldata <- data.frame(condition = donor_condition)       # one condition label per donor, aligned to pseudobulk's columns
dds <- DESeqDataSetFromMatrix(pseudobulk, coldata, ~ condition)
dds <- DESeq(dds)
disp_vec <- dispersions(dds); mu_vec <- rowMeans(counts(dds, normalized = TRUE))
keep <- is.finite(disp_vec) & is.finite(mu_vec) & mu_vec > 0
disp_vec <- disp_vec[keep]; mu_vec <- mu_vec[keep]

res <- safe_ssize(function() ssizeRNA_vary(nGenes = length(mu_vec), pi0 = 0.95,
                                            mu = mu_vec, disp = disp_vec, fc = fc, fdr = fdr,
                                            power = power, maxN = maxN))
n <- require_reachable_n(res, maxN)                         # minimum DONORS per group -- NOT cells

stopifnot(!anyNA(donor_condition), ncol(pseudobulk) == length(cell_counts),
          length(mu_vec) > 0, is.finite(median(disp_vec)), n >= 2, n <= maxN)
cat(sprintf("OK: pseudobulk median dispersion %.3f; minimum donors per group = %d (fc=%.2f, achieved power %.3f)\n",
            median(disp_vec), n, fc, res$ssize[, "power"]))
