# Multiple testing: FDR vs FWER, dependence (BH vs BY), q-value/pi0, local FDR, IHW
# Reference: R stats (base), qvalue 2.34+, IHW 1.30+ | Verify API if version differs
#
# Demonstrates the error-rate choices that matter in genomics discovery and the levers
# (pi0, covariate weighting) that buy back power. Python statsmodels note: the DEFAULT
# method is 'hs' (Holm-Sidak), NOT BH -- always pass method= explicitly.

set.seed(20260528)
n_genes <- 10000; n_de <- 500
pvalues <- c(rbeta(n_de, 0.3, 5), runif(n_genes - n_de))   # 5% true DE (small p), rest null
is_de  <- c(rep(TRUE, n_de), rep(FALSE, n_genes - n_de))

# ---------------------------------------------------------------------------
# 1. FWER (Bonferroni/Holm) vs FDR (BH) vs FDR-under-dependence (BY), averaged over draws
# ---------------------------------------------------------------------------
# The realized FDP of ONE draw is noisy (BH's FDP varies a lot run to run), so a single seed can show
# BH above nominal by chance. FDR is an expectation: average the table over replicate draws.
methods <- c(none = NA, bonferroni = 'bonferroni', holm = 'holm', BH = 'BH', BY = 'BY')
n_rep <- 50
acc <- array(0, dim = c(length(methods), 3), dimnames = list(names(methods), c('significant', 'true_pos', 'false_pos')))
fdp <- matrix(NA_real_, n_rep, length(methods), dimnames = list(NULL, names(methods)))
for (r in seq_len(n_rep)) {
  p_r <- c(rbeta(n_de, 0.3, 5), runif(n_genes - n_de))
  for (i in seq_along(methods)) {
    padj <- if (is.na(methods[i])) p_r else p.adjust(p_r, method = methods[i])
    sig <- padj < 0.05
    acc[i, ] <- acc[i, ] + c(sum(sig), sum(sig & is_de), sum(sig & !is_de))
    fdp[r, i] <- sum(sig & !is_de) / max(sum(sig), 1)
  }
}
tab <- data.frame(method = names(methods), acc / n_rep, mean_fdp = round(colMeans(fdp), 3), row.names = NULL)
print(tab, row.names = FALSE, digits = 4)   # per-replicate means over n_rep draws; mean_fdp is the realized FDR
# BH is the discovery default (valid under independence/PRDS); BY is conservative but valid
# under arbitrary/negative dependence.

# ---------------------------------------------------------------------------
# 2. q-value: estimate pi0 (true-null proportion) for more power; local FDR per feature
# ---------------------------------------------------------------------------
# Do not load qvalue in this process: requireNamespace('qvalue') can crash this Windows R
# session during teardown. The isolated helper records qvalue's output before its worker exits.
source('scripts/qvalue_safe.R')
qres <- qvalue_safe(pvalues)
if (identical(qres$method, 'qvalue')) {
  cat(sprintf('\nq-value: pi0 = %.3f, discoveries at q<0.05 = %d\n',
              qres$pi0, qres$discoveries))
} else {
  cat(sprintf('\nq-value did not produce output; %s = %d discoveries\n',
              qres$method, qres$discoveries))
}

# ---------------------------------------------------------------------------
# 3. IHW: weight hypotheses by an independent covariate (must be null-independent)
# ---------------------------------------------------------------------------
if (requireNamespace('IHW', quietly = TRUE)) {
  # The IHW solver can SIGSEGV in a child R process. Reuse the tested wrapper rather than
  # maintaining a second process-management implementation in this example.
  source('scripts/ihw_safe.R')
  mean_expr <- rgamma(n_genes, shape = 2, rate = 0.5)        # covariate independent of null p
  ihw_res <- ihw_safe(pvalues, mean_expr, alpha = 0.05)
  bh_discoveries <- sum(p.adjust(pvalues, 'BH') < 0.05)
  if (identical(ihw_res$method, 'IHW')) {
    cat(sprintf('IHW discoveries at FDR 0.05 = %d (vs BH = %d; attempt %d)\n',
                sum(ihw_res$padj < 0.05), bh_discoveries, ihw_res$attempts))
  } else {
    cat(sprintf('IHW did not complete after %d attempts (solver crash); falling back to BH = %d discoveries\n',
                ihw_res$attempts, bh_discoveries))
  }
}

# ---------------------------------------------------------------------------
# 4. Python equivalent (mind the default method)
# ---------------------------------------------------------------------------
# from statsmodels.stats.multitest import multipletests
# rej, padj, _, _ = multipletests(pvalues, alpha=0.05, method='fdr_bh')   # BH -- NOT the default
# rej, padj, _, _ = multipletests(pvalues, alpha=0.05, method='fdr_by')   # BY under dependence

# GWAS: genome-wide significance ~5e-8 (Dudbridge & Gusnanto 2008 derived ~7.2e-8);
# the GWAS test machinery lives in population-genetics/association-testing.
