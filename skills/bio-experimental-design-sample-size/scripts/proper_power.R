# Purpose: simulate power by replicate number from a pilot count matrix with PROPER
#          (estParam -> RNAseq.SimOptions.2grp -> runSims -> comparePower).
# Inputs:  pilot counts CSV (genes x samples, gene IDs in column 1), then optional
#          positional args: reps (comma list, default 3,6,10,20), nsims (20), fc (1.5),
#          max_genes (0 = all genes), seed (20260918).
# Usage:   r.sh proper_power.R pilot_counts.csv [3,6,10,20] [20] [1.5] [0] [20260918]
# Checked: PROPER 1.38.0, R 4.4.3 / Bioconductor 3.20.
suppressPackageStartupMessages(library(PROPER))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: proper_power.R pilot_counts.csv [reps] [nsims] [fc] [max_genes] [seed]")
arg <- function(i, default) if (length(args) >= i) args[i] else default
counts_file <- args[1]
reps      <- as.integer(strsplit(arg(2, "3,6,10,20"), ",")[[1]])
nsims     <- as.integer(arg(3, 20))
fc        <- as.numeric(arg(4, 1.5))
max_genes <- as.integer(arg(5, 0))
seed      <- as.integer(arg(6, 20260918))
set.seed(seed)

pilot_counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
if (max_genes > 0) pilot_counts <- pilot_counts[seq_len(min(max_genes, nrow(pilot_counts))), , drop = FALSE]
counts_mat <- as.matrix(pilot_counts)
oldClass(counts_mat) <- "matrix"        # work around PROPER::estParam's `class(X) %in% c(...)` check,
                                         # which errors on R >= 4.0's matrix/array dual class (verified PROPER 1.38.0)
params <- estParam(counts_mat, type = 1)
sim.opts <- RNAseq.SimOptions.2grp(ngenes = nrow(counts_mat), seqDepth = params$seqDepth,
                                    lBaselineExpr = params$lmean, lOD = params$lOD,
                                    p.DE = 0.05, lfc = log2(fc), sim.seed = seed)
simres <- runSims(Nreps = reps, nsims = nsims, sim.opts = sim.opts, DEmethod = "DESeq2")
# comparePower counts a DE gene as a target only if abs(lfc) > delta (STRICT). runSims planted every DE gene at
# exactly lfc = log2(fc), so delta = log2(fc) leaves ZERO target genes and power.marginal is all NaN.
# Put delta just below the planted fold change:
powres <- comparePower(simres, alpha.type = "fdr", alpha.nominal = 0.05,
                       stratify.by = "expr", target.by = "lfc", delta = log2(fc) - 0.01)
names(powres)                           # 16-field list; there is no `powerAveraged` (`$` on a missing name returns NULL silently)
powres$Nreps1                           # replicate numbers, same order as the rows below
marginal <- rowMeans(powres$power.marginal, na.rm = TRUE)   # power.marginal is Nreps x nsims; marginal power per Nreps
print(marginal)
summaryPower(powres)                    # PROPER's own table: nominal vs actual FDR, marginal power, avg TD/FD per Nreps

stopifnot(length(marginal) == length(reps), !is.null(powres$Nreps1),
          any(is.finite(marginal)), all(marginal >= 0 & marginal <= 1, na.rm = TRUE))
cat("OK: marginal power per Nreps =", paste(sprintf("%d:%.3f", reps, marginal), collapse = "  "), "\n")
