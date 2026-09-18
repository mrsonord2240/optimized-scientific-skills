# Sample size for genomics: FDR-aware NB sizing, pilot dispersions, verifying a fixed n
# Reference: checked on ssizeRNA 1.3.3, DESeq2 1.46.0 (R 4.4.3 / Bioconductor 3.20), 2026-09-18
#
# Demonstrates: ssizeRNA_single (note: `m` is the pseudo sample size for SIMULATION,
# default 200 -- NOT the number of DE genes or n per group; `res$ssize` is a 1x3 matrix
# (pi0, ssize, power) -- always index the column you want), estimating dispersion from a
# pilot, check.power() to confirm a budget-fixed n meets the target FDR, and what to do
# when no n within maxN reaches the target.

suppressPackageStartupMessages(library(ssizeRNA))
set.seed(20260528)

# ---------------------------------------------------------------------------
# 1. FDR-aware sample size (single mean/dispersion for all genes)
# ---------------------------------------------------------------------------
# nGenes: total genes; pi0: proportion NON-DE; m: pseudo sample size for simulation (default 200);
# mu: mean NORMALIZED count -- use a realistic value (order 100-500 for typical bulk RNA-seq
#     depth), not a small toy number, or the search can run past maxN and come back NA.
res <- ssizeRNA_single(nGenes = 20000, pi0 = 0.95, m = 200,
                       mu = 200, disp = 0.2, fc = 1.5, fdr = 0.05, power = 0.80, maxN = 200)
n <- res$ssize[, "ssize"]
if (is.na(n)) stop("no n <= maxN reaches the target; raise maxN or revise fc/dispersion")
cat(sprintf('Minimum n per group (1.5-fold, mu=200, disp=0.2, FDR 0.05, 80%% power): %d (achieved power %.3f)\n',
            n, res$ssize[, "power"]))

# Sensitivity to fold change (the dominant lever)
for (fc in c(1.5, 2, 3)) {
  r <- ssizeRNA_single(nGenes = 20000, pi0 = 0.95, m = 200, mu = 200, disp = 0.2,
                       fc = fc, fdr = 0.05, power = 0.80, maxN = 200)
  n_fc <- r$ssize[, "ssize"]
  cat(sprintf('fc=%.1f -> n=%s per group\n', fc, if (is.na(n_fc)) "NA (raise maxN)" else n_fc))
}

# ---------------------------------------------------------------------------
# 1b. Reconcile the ">=6 replicates" rule of thumb against this calculation
# ---------------------------------------------------------------------------
# Schurch 2016's ">=6" is an empirical RECOVERY benchmark averaged over a real spectrum of
# fold changes (most genes change by less than any single target FC). The ssizeRNA number
# above answers a stricter, different question: 80% MARGINAL power at one FIXED minimum
# fold change. They are not interchangeable -- print both so neither is mistaken for the other.
cat(sprintf('\nSchurch 2016 empirical floor: >=6 replicates recovers most true DE across a realistic FC spectrum.\n'))
cat(sprintf('This calculation (80%% power at a FIXED %.1f-fold, FDR 0.05): n=%d per group.\n', 1.5, n))

# ---------------------------------------------------------------------------
# 2. Estimate dispersion from a pilot (defensible input) and size with vectors
# ---------------------------------------------------------------------------
if (requireNamespace('DESeq2', quietly = TRUE)) {
  suppressPackageStartupMessages(library(DESeq2))
  # dds  <- DESeqDataSetFromMatrix(pilot_counts, pilot_coldata, ~ condition)
  # dds  <- DESeq(dds)
  # disp_vec <- dispersions(dds); mu_vec <- rowMeans(counts(dds, normalized = TRUE))
  # ssizeRNA_vary(nGenes = length(disp_vec), pi0 = 0.95, mu = mu_vec, disp = disp_vec,
  #               fc = 1.5, fdr = 0.05, power = 0.80, maxN = 200)$ssize[, "ssize"]
  # ssizeRNA_vary needs per-gene VECTORS -- passing scalars raises
  # 'integrate(): non-finite function value' regardless of the values (use ssizeRNA_single instead; see above).
  cat('\nUse DESeq2::dispersions(dds) + normalized means from a pilot as VECTORS to ssizeRNA_vary.\n')
}

# ---------------------------------------------------------------------------
# 3. Verify a budget-fixed n: average power and TRUE realized FDR (here m = n per group)
# ---------------------------------------------------------------------------
cp <- check.power(nGenes = 20000, pi0 = 0.95, m = 6, mu = 200, disp = 0.2,
                  fc = 1.5, fdr = 0.05, sims = 50)
if (is.nan(cp$fdr_bh_ave)) {
  cat(sprintf('\nAt n=6/group: BH average power = %.3f, true FDR = NaN (zero discoveries, NOT "FDR unknown").\n',
              cp$pow_bh_ave))
} else {
  cat(sprintf('\nAt n=6/group: BH average power = %.3f, true FDR = %.3f\n', cp$pow_bh_ave, cp$fdr_bh_ave))
}

# ---------------------------------------------------------------------------
# 4. Assay floors (NOT targets) -- defensible only after pilot/literature dispersion supports them
# ---------------------------------------------------------------------------
floors <- data.frame(
  assay = c('Bulk RNA-seq', 'scRNA-seq (population DE)', 'ATAC-seq', 'ChIP-seq',
            'Proteomics (DIA/TMT)', 'Methylation (WGBS)'),
  min_replicates = c(3, 3, 2, 2, 3, 4),                 # floors under low dispersion + large effects
  for_small_effects = c('6-12', '6+ donors', '4-6', '3-4', '6-10', '8-12'),
  note = c('Schurch 2016: >=6 recovers most true DE across a spectrum -- NOT the fixed-FC n above (see 1b)',
           'donors, not cells, set power (Squair 2021)',
           'library complexity floor', 'IDR reproducibility (ENCODE)',
           'floor assumes single-marker validation, NOT proteome-wide FDR control -- adjust alpha for the panel',
           'high per-CpG variance'))
print(floors, row.names = FALSE)
cat('\n# Count biological replicates, not measurements; add 10-20% for sample failures\n')
cat('# (for human-donor cohorts, the failure margin must match the approved protocol).\n')
