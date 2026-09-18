# Reference: coloc 5.2.3+, susieR 0.12.35+ | Verify API if version differs
## Multi-causal coloc.susie pipeline with LD-z-score consistency diagnostics
##
## Demonstrates:
##   1. Self-consistent genotype -> summary-stats -> LD generation. coloc.susie
##      routes through SuSiE fine-mapping, so its LD matrix and z-scores MUST
##      come from the same genotype process -- a hand-built toy LD matrix
##      generated independently of the simulated betas/SEs is not consistent
##      with any real genotype structure, fails estimate_s_rss() below, and
##      makes runsusie() abort with "the estimated prior variance is
##      unreasonably large". This example instead simulates individual-level
##      genotypes once and derives both the summary statistics and the LD
##      matrix from that same matrix (SKILL.md's "Use in-sample LD when at
##      all possible" guidance).
##   2. estimate_s_rss diagnostic (LD reference must match z-scores)
##   3. SuSiE per-trait credible sets with L = 10 max signals
##   4. coloc.susie per (CS1, CS2) pair
##   5. Reporting per-pair PP.H4 with credible-set lead SNPs
##
## GWAS has two independent causal signals (allelic heterogeneity); only one
## is shared with the eQTL.

library(coloc)
library(susieR)

set.seed(303)
n_ind <- 4000
n_snps <- 300
positions <- sort(sample(30000000:31000000, n_snps))

# --- Simulate genotype dosages with AR(1)-like LD ---
rho <- 0.85
Sigma <- rho^abs(outer(1:n_snps, 1:n_snps, '-'))
L <- chol(Sigma)
G <- matrix(rnorm(n_ind * n_snps), n_ind, n_snps) %*% L
G <- scale(G)   # standardised dosage-like genotype

## Two independent causal SNPs in GWAS; only causal1 is shared with eQTL
causal1 <- which.min(abs(positions - 30300000))
causal2 <- which.min(abs(positions - 30700000))
cat('PLANTED TRUTH: causal1 idx', causal1, '(shared GWAS+eQTL), causal2 idx', causal2, '(GWAS-only)\n')

y_gwas <- 0.25 * G[, causal1] + 0.20 * G[, causal2] + rnorm(n_ind, 0, 1)
y_eqtl <- 0.40 * G[, causal1] + rnorm(n_ind, 0, 1)   # only causal1 shared with eQTL

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

## In-sample LD -- derived from the SAME genotype matrix used for both traits.
## coloc::runsusie matches the `snp` vector to dimnames(LD); unnamed LD errors out.
ld_matrix <- cor(G)
dimnames(ld_matrix) <- list(snp_ids, snp_ids)

gwas_N <- n_ind
eqtl_N <- n_ind

## Critical diagnostic: lambda > 0.05 means LD does not match z-scores -> abort
z_gwas <- gwas_ss$beta / gwas_ss$se
lam_gwas <- susieR::estimate_s_rss(z = z_gwas, R = ld_matrix, n = gwas_N)
cat(sprintf('GWAS estimate_s_rss lambda = %.4f\n', lam_gwas))
if (lam_gwas > 0.05) stop('GWAS LD-z mismatch; abort or use coloc.abf')

z_eqtl <- eqtl_ss$beta / eqtl_ss$se
lam_eqtl <- susieR::estimate_s_rss(z = z_eqtl, R = ld_matrix, n = eqtl_N)
cat(sprintf('eQTL estimate_s_rss lambda = %.4f\n', lam_eqtl))
if (lam_eqtl > 0.05) stop('eQTL LD-z mismatch; abort or use coloc.abf')

gwas_data <- list(beta = gwas_ss$beta, varbeta = gwas_ss$se^2, snp = snp_ids, position = positions,
                   type = 'cc', s = 0.3, N = gwas_N, LD = ld_matrix)
eqtl_data <- list(beta = eqtl_ss$beta, varbeta = eqtl_ss$se^2, snp = snp_ids, position = positions,
                   type = 'quant', sdY = sd(y_eqtl), N = eqtl_N, LD = ld_matrix)

## L = 10 = max number of credible sets to detect; SuSiE returns fewer when fewer are supported
s_gwas <- runsusie(gwas_data, L = 10)
s_eqtl <- runsusie(eqtl_data, L = 10)

n_cs_gwas <- if (!is.null(summary(s_gwas)$cs)) nrow(summary(s_gwas)$cs) else 0
n_cs_eqtl <- if (!is.null(summary(s_eqtl)$cs)) nrow(summary(s_eqtl)$cs) else 0
cat(sprintf('\nGWAS credible sets: %d (planted truth: 2) | eQTL credible sets: %d (planted truth: 1)\n',
            n_cs_gwas, n_cs_eqtl))

res_susie <- coloc.susie(s_gwas, s_eqtl)

if (is.null(res_susie$summary) || nrow(res_susie$summary) == 0) {
    cat('\nNo overlapping credible sets between traits -> no colocalization signal\n')
} else {
    cat('\nPer-(CS1, CS2) colocalization PP:\n')
    print(res_susie$summary[, c('idx1', 'idx2', 'hit1', 'hit2',
                                  'PP.H3.abf', 'PP.H4.abf')])
    best <- res_susie$summary[which.max(res_susie$summary$PP.H4.abf), ]
    cat(sprintf('\nBest CS pair: GWAS CS%d (%s) x eQTL CS%d (%s) | PP.H4 = %.3f (planted shared causal: rs%d)\n',
                 best$idx1, best$hit1, best$idx2, best$hit2, best$PP.H4.abf, causal1))
}
