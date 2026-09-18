# Reference: coloc 5.2.3+, susieR 0.12.35+ | Verify API if version differs
## SuSiE-coloc basic case: one shared causal variant
##
## coloc.susie routes through SuSiE fine-mapping, so its LD matrix and its
## z-scores must come from the SAME genotype process -- an LD matrix built by
## hand (e.g. an exponential-decay toy) is not automatically consistent with
## betas/SEs simulated independently, and susieR's own estimate_s_rss()
## diagnostic (see coloc_susie_multicausal.R) will correctly reject it. This
## example instead simulates individual-level genotypes once and derives both
## the summary statistics and the LD matrix from that same matrix, exactly as
## SKILL.md's "Use in-sample LD when at all possible" guidance recommends.

library(coloc)
library(susieR)

set.seed(42)
n_ind <- 4000
n_snps <- 300
positions <- sort(sample(30000000:31000000, n_snps))

# --- Simulate genotype dosages with AR(1)-like LD ---
rho <- 0.85
Sigma <- rho^abs(outer(1:n_snps, 1:n_snps, '-'))
L <- chol(Sigma)
G <- matrix(rnorm(n_ind * n_snps), n_ind, n_snps) %*% L
G <- scale(G)   # standardised dosage-like genotype

causal <- which.min(abs(positions - 30500000))
cat('PLANTED TRUTH: shared causal SNP index', causal, '\n')

# GWAS and eQTL phenotypes both driven by the same causal SNP
y_gwas <- 0.30 * G[, causal] + rnorm(n_ind, 0, 1)
y_eqtl <- 0.40 * G[, causal] + rnorm(n_ind, 0, 1)

# Per-SNP marginal regression (as a real GWAS/eQTL summary-stats pipeline would produce)
get_sumstats <- function(y, G) {
  beta <- se <- numeric(ncol(G))
  for (j in 1:ncol(G)) {
    fit <- summary(lm(y ~ G[, j]))$coefficients
    beta[j] <- fit[2, 1]
    se[j] <- fit[2, 2]
  }
  list(beta = beta, se = se)
}
gwas_ss <- get_sumstats(y_gwas, G)
eqtl_ss <- get_sumstats(y_eqtl, G)

snp_ids <- paste0('rs', 1:n_snps)

# In-sample LD -- derived from the SAME genotype matrix used for both traits.
# coloc::runsusie matches the `snp` vector to dimnames(LD); unnamed LD errors out.
ld_matrix <- cor(G)
dimnames(ld_matrix) <- list(snp_ids, snp_ids)

# --- Format for SuSiE ---
gwas_data <- list(
  beta = gwas_ss$beta, varbeta = gwas_ss$se^2,
  snp = snp_ids, position = positions,
  type = 'cc', s = 0.3, N = n_ind, LD = ld_matrix
)
eqtl_data <- list(
  beta = eqtl_ss$beta, varbeta = eqtl_ss$se^2,
  snp = snp_ids, position = positions,
  type = 'quant', sdY = sd(y_eqtl), N = n_ind, LD = ld_matrix
)

# --- Run SuSiE on each dataset ---
# L: Maximum number of causal variants to consider
# L = 10 is a reasonable default even for a single true signal; SuSiE
# returns fewer credible sets when fewer are supported by the data.
susie_gwas <- runsusie(gwas_data, L = 10)
susie_eqtl <- runsusie(eqtl_data, L = 10)

n_gwas_cs <- if (is.null(summary(susie_gwas)$cs)) 0 else nrow(summary(susie_gwas)$cs)
n_eqtl_cs <- if (is.null(summary(susie_eqtl)$cs)) 0 else nrow(summary(susie_eqtl)$cs)
cat('GWAS credible sets found:', n_gwas_cs, '(planted truth: 1)\n')
cat('eQTL credible sets found:', n_eqtl_cs, '(planted truth: 1)\n')

# --- Run SuSiE-coloc ---
result <- coloc.susie(susie_gwas, susie_eqtl)

if (!is.null(result$summary)) {
  cat('\nSuSiE-coloc results:\n')
  print(result$summary[, c('hit1', 'hit2', 'PP.H3.abf', 'PP.H4.abf')])
  # hit1, hit2: lead SNP of the GWAS / eQTL credible set in this pair
  # PP.H4.abf: Posterior probability of shared causal variant for this pair
  best <- result$summary[which.max(result$summary$PP.H4.abf), ]
  cat(sprintf('\nBest CS pair: %s x %s | PP.H4 = %.3f (planted causal: rs%d)\n',
              best$hit1, best$hit2, best$PP.H4.abf, causal))
} else {
  cat('\nNo overlapping credible sets found\n')
}
