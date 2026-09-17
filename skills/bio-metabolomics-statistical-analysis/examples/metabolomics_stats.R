# Reference: ropls 1.34+ | Verify API if version differs
# Permutation-validated OPLS-DA on synthetic metabolomics data, demonstrating that
# (1) only pQ2 -- not R2Y or the score plot -- licenses a discriminant claim, and
# (2) the scaling choice (Pareto vs unit-variance) is a hypothesis that changes the VIP list.
library(ropls)

set.seed(1)
n_per_group <- 20
n_features <- 300
n_true <- 15                                  # only the first 15 features actually differ

# Synthetic data is already homoscedastic and roughly normal; real MS intensities are
# right-skewed and heteroscedastic, so log/glog-transform them BEFORE opls() in practice.
group <- factor(rep(c('control', 'case'), each = n_per_group))
intensities <- matrix(rnorm(2 * n_per_group * n_features, mean = 10, sd = 1),
                      nrow = 2 * n_per_group, ncol = n_features)
colnames(intensities) <- paste0('M', seq_len(n_features))
case_rows <- group == 'case'
intensities[case_rows, seq_len(n_true)] <- intensities[case_rows, seq_len(n_true)] + 1.2

# A NULL model: same data, labels permuted. OPLS-DA will still separate it in p>>n,
# so a clean score plot here is the geometry, not signal -- pQ2 is the honest check.
null_group <- factor(sample(as.character(group)))

# ropls's own cross-validated significance test on the first predictive component can
# reject it and silently return a 0-row summaryDF / empty model -- with info.txtC='none'
# this produces NO warning or error, yet class(model) still reads "opls" and getVipVn()
# returns a length-0 vector instead of erroring (measured at ~40% of runs at this n/p).
# ALWAYS check nrow(getSummaryDF()) before trusting a fit; never read R2/Q2/VIP off an
# unchecked model. See SKILL.md's "OPLS-DA silently returns an empty model" entry.
run_model <- function(y, scaleC) {
    m <- opls(intensities, y, predI = 1, orthoI = NA, scaleC = scaleC,
              permI = 1000, crossvalI = 7, fig.pdfC = 'none', info.txtC = 'none')
    if (nrow(getSummaryDF(m)) > 0) return(list(model = m, type = 'OPLS-DA'))
    # Empty OPLS-DA: fall back to PLS-DA (orthoI=0) -- identical predictive power, no
    # orthogonal-component significance gate to fail; verified 0/30 failures on both
    # signal-bearing and pure-noise synthetic data at this n/p (see the Skill's fix log).
    m2 <- opls(intensities, y, predI = 1, orthoI = 0, scaleC = scaleC,
               permI = 1000, crossvalI = 7, fig.pdfC = 'none', info.txtC = 'none')
    if (nrow(getSummaryDF(m2)) > 0) return(list(model = m2, type = 'PLS-DA (OPLS-DA fallback)'))
    stop('Neither OPLS-DA nor the PLS-DA fallback produced a usable model: the first ',
         'predictive component was not significant under ropls\' own cross-validated ',
         'criterion. Report that no multivariate separation was detected -- do not force a model.')
}

cat('=== Real labels, Pareto scaling ===\n')
fit_real_pareto <- run_model(group, 'pareto')
real_pareto <- fit_real_pareto$model
cat('Model type actually fit:', fit_real_pareto$type, '\n')
print(getSummaryDF(real_pareto)[, c('R2X(cum)', 'R2Y(cum)', 'Q2(cum)', 'pR2Y', 'pQ2')])

cat('\n=== Permuted (null) labels, Pareto scaling -- expect high R2Y, failed pQ2 ===\n')
fit_null_pareto <- run_model(null_group, 'pareto')
null_pareto <- fit_null_pareto$model
cat('Model type actually fit:', fit_null_pareto$type, '\n')
print(getSummaryDF(null_pareto)[, c('R2X(cum)', 'R2Y(cum)', 'Q2(cum)', 'pR2Y', 'pQ2')])

cat('\n=== Real labels, unit-variance scaling -- compare the VIP ranking ===\n')
fit_real_uv <- run_model(group, 'standard')
real_uv <- fit_real_uv$model
cat('Model type actually fit:', fit_real_uv$type, '\n')
print(getSummaryDF(real_uv)[, c('R2X(cum)', 'R2Y(cum)', 'Q2(cum)', 'pR2Y', 'pQ2')])

vip_pareto <- getVipVn(real_pareto)
vip_uv <- getVipVn(real_uv)
top10_pareto <- names(sort(vip_pareto, decreasing = TRUE))[1:10]
top10_uv <- names(sort(vip_uv, decreasing = TRUE))[1:10]
overlap <- length(intersect(top10_pareto, top10_uv))
cat('\nTop-10 VIP overlap between Pareto and UV scaling:', overlap, 'of 10\n')
cat('Scaling-fragile result if the two lists diverge.\n')

# pQ2 (the permutation p-value for Q2) is the licensing gate -- it, not R2Y or the score
# plot, exposes a model built on noise. Q2's magnitude is the effect size and should also
# clear the Triba >0.5 heuristic to be worth reporting.
summ <- getSummaryDF(real_pareto)
licensed <- summ[['pQ2']] < 0.05
cat(sprintf('\nReal model: Q2(cum)=%.2f, pQ2=%.3f -> permutation-licensed: %s\n',
            summ[['Q2(cum)']], summ[['pQ2']], licensed))
