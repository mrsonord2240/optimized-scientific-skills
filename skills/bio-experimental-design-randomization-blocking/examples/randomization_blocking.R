# Randomization, blocking, and the experimental unit
# Reference: designit 0.5+, lme4 1.1-35+, lmerTest 3.1+ | Verify API if version differs
#
# Demonstrates: identifying the experimental unit and aggregating observational units,
# restricted (block) randomization with run-order randomization, and the mixed-model
# random-effects structure that matches a nested / split-plot design.

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
})

# ---------------------------------------------------------------------------
# 1. The experimental unit defines n: aggregate cells (observational units) to donors
# ---------------------------------------------------------------------------
set.seed(20260528)
n_donors_per_group <- 4
cells_per_donor <- 200

cells <- expand.grid(cell = seq_len(cells_per_donor),
                     donor = paste0('D', seq_len(2 * n_donors_per_group)))
cells$condition <- ifelse(as.integer(sub('D', '', cells$donor)) <= n_donors_per_group,
                         'ctrl', 'treat')
# donor-level mean shift + cell-level noise
donor_effect <- setNames(rnorm(length(unique(cells$donor)), 0, 0.5), unique(cells$donor))
cells$measurement <- donor_effect[cells$donor] +
  ifelse(cells$condition == 'treat', 0.3, 0) + rnorm(nrow(cells), 0, 1)

# Correct unit of inference: one value per donor (the experimental unit), not per cell.
eu_level <- cells |>
  group_by(donor, condition) |>
  summarise(value = mean(measurement), .groups = 'drop')
cat('Experimental units (n) per group:',
    paste(table(eu_level$condition), collapse = ' / '), '\n')   # 4 / 4, NOT 800 / 800

# Test on experimental-unit values (n = 8 donors), not on 1600 pseudoreplicated cells.
t.test(value ~ condition, data = eu_level)

# ---------------------------------------------------------------------------
# 2. Restricted (block) randomization + run-order randomization
# ---------------------------------------------------------------------------
# 24 samples processed 8/day over 3 days; randomize condition WITHIN day (block),
# and randomize the processing order so position is not confounded with condition.
units <- data.frame(id = sprintf('S%02d', 1:24),
                    day = rep(c('day1', 'day2', 'day3'), each = 8))
units$condition <- ave(units$id, units$day, FUN = function(ids)
  sample(rep(c('ctrl', 'treat'), length.out = length(ids))))
units$run_order <- sample(nrow(units))

# Confirm balance by assertion, not by eye: each day must hold equal ctrl/treat
# (day orthogonal to condition), and run_order must be a valid within-day permutation.
bal <- table(units$day, units$condition)
print(bal)
stopifnot(
  "each day must have exactly 4 ctrl and 4 treat" = all(bal == 4),
  "run_order must be a permutation of 1:24" = setequal(units$run_order, seq_len(nrow(units)))
)
cat('Balance assertions passed.\n')

# ---------------------------------------------------------------------------
# 3. Mixed model matching a nested design (cells within donor)
# ---------------------------------------------------------------------------
# Random intercept for donor accounts for the within-donor correlation that makes
# cells pseudoreplicates; condition stays a fixed effect.
fit_nested <- lmer(measurement ~ condition + (1 | donor), data = cells)
print(anova(fit_nested))   # Satterthwaite df via lmerTest; effective n driven by donors

# Split-plot pattern (whole plot = processing day/run, sub-plot = condition):
#   lmer(response ~ condition + (1 | day/sample), data = df)
# The whole-plot factor must be tested against whole-plot error, never sub-plot error.

# ---------------------------------------------------------------------------
# 4. Split-plot Type-I error demonstration: a flat model IS anti-conservative for a
#    whole-plot fixed effect; the split-plot mixed model is NOT.
# ---------------------------------------------------------------------------
# temp is a WHOLE-PLOT factor (assigned per run; true effect = 0, i.e. the null).
# genotype is a SUB-PLOT factor (randomized within run; true effect = 0.6).
# Repeat the simulated experiment many times and count how often each model's temp
# p-value falls below 0.05 under the null -- the nominal Type-I error rate is 5%.
n_runs      <- 6
n_sub       <- 4
sigma_run   <- 1.5   # between-run (whole-plot) SD -- the nuisance a flat model ignores
sigma_resid <- 0.5
true_genotype_effect <- 0.6
true_temp_effect     <- 0    # null: no real whole-plot effect
nsim  <- 400
alpha <- 0.05

one_sim <- function() {
  runs <- data.frame(run = factor(1:n_runs), temp = rep(c('low', 'high'), each = n_runs / 2))
  runs$temp <- sample(runs$temp)
  run_effect <- setNames(rnorm(n_runs, 0, sigma_run), runs$run)
  df <- do.call(rbind, lapply(seq_len(n_runs), function(r) {
    genotype <- sample(rep(c('WT', 'KO'), each = n_sub / 2))
    data.frame(run = runs$run[r], temp = runs$temp[r], genotype = genotype)
  }))
  df$y <- true_temp_effect * (df$temp == 'high') + true_genotype_effect * (df$genotype == 'KO') +
    run_effect[df$run] + rnorm(nrow(df), 0, sigma_resid)

  flat    <- lm(y ~ temp + genotype, data = df)                 # WRONG: ignores run clustering
  correct <- lmer(y ~ temp + genotype + (1 | run), data = df)   # CORRECT: whole-plot error stratum
  c(flat = coef(summary(flat))['templow', 'Pr(>|t|)'],
    correct = coef(summary(correct))['templow', 'Pr(>|t|)'])
}

set.seed(20260917)
res <- t(replicate(nsim, one_sim()))
rate_flat    <- mean(res[, 'flat']    < alpha)
rate_correct <- mean(res[, 'correct'] < alpha)
cat(sprintf('Flat lm() rejection rate for the null whole-plot effect:    %.3f  (anti-conservative if >> %.2f)\n',
            rate_flat, alpha))
cat(sprintf('Split-plot lmer(...+(1|run)) rejection rate:                %.3f  (should be close to %.2f)\n',
            rate_correct, alpha))
stopifnot(
  "flat model must be clearly anti-conservative (>3x nominal alpha)" = rate_flat > 3 * alpha,
  "correct split-plot model must be within a reasonable band of nominal alpha" =
    rate_correct > alpha / 2.5 && rate_correct < alpha * 2.5,
  "flat model's Type-I error must exceed the correct model's by a wide margin" =
    rate_flat > rate_correct + 0.10
)
