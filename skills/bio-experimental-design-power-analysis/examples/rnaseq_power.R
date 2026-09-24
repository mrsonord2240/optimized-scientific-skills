# Power analysis for RNA-seq: closed-form vs simulation, per-gene power
# Reference: RNASeqPower 1.42+, PROPER 1.34+ | Verify API if version differs
#
# Demonstrates: closed-form NB power (both directions), why a single CV is only an
# approximation, simulation-based marginal power at a target FDR, and the depth-vs-replicate
# tradeoff. Observed/post-hoc power is deliberately NOT computed -- it is uninformative.

suppressPackageStartupMessages(library(RNASeqPower))

# Guard RNASeqPower's permissive argument handling before every prospective call.
validate_rnapower_inputs <- function(depth, cv, effect = NULL, alpha = 0.05,
                                     n = NULL, power = NULL) {
  check_scalar <- function(x, name, predicate, rule) {
    if (length(x) != 1L || !is.finite(x) || !predicate(x)) stop(sprintf("%s must be one finite scalar %s", name, rule), call. = FALSE)
  }
  check_scalar(depth, "depth", function(x) x > 0, "> 0")
  check_scalar(cv, "cv", function(x) x > 0, "> 0")
  check_scalar(alpha, "alpha", function(x) x > 0 && x < 1, "in (0, 1)")
  if (!is.null(n)) check_scalar(n, "n", function(x) x > 0, "> 0")
  if (!is.null(power)) check_scalar(power, "power", function(x) x > 0 && x < 1, "in (0, 1)")
  if (!is.null(effect)) check_scalar(effect, "effect", function(x) x > 0 && x != 1, "> 0 and != 1")
  if (sum(vapply(list(n, power, effect), is.null, logical(1))) != 1L) stop("Provide exactly two of n, power, and effect; RNASeqPower solves the omitted quantity.", call. = FALSE)
  invisible(TRUE)
}
checked_rnapower <- function(depth, cv, effect = NULL, alpha = 0.05, n = NULL, power = NULL) {
  validate_rnapower_inputs(depth, cv, effect, alpha, n, power)
  do.call(RNASeqPower::rnapower, Filter(Negate(is.null), list(depth = depth, n = n, cv = cv, effect = effect, alpha = alpha, power = power)))
}

assess_realized_fdr <- function(summary_table, target_fdr = 0.05, tolerance = 0) {
  actual <- as.numeric(summary_table[, "Actual FDR"])
  data.frame(replicates_per_group = summary_table[, "SS1"], actual_fdr = actual,
             marginal_power = as.numeric(summary_table[, "Marginal power"]),
             accepted = is.finite(actual) & actual <= target_fdr + tolerance,
             decision = ifelse(is.finite(actual) & actual <= target_fdr + tolerance, "ACCEPT", "REJECT"))
}

# ---------------------------------------------------------------------------
# 1. Closed-form NB power -- solves for whichever of n / power is omitted
# ---------------------------------------------------------------------------
# depth = reads/gene; cv = biological coefficient of variation; effect = fold change
cat('Power, n=3, 2-fold, CV=0.4:', round(checked_rnapower(depth = 20, n = 3, cv = 0.4, effect = 2, alpha = 0.05), 3), '\n')
cat('Power, n=6, 2-fold, CV=0.4:', round(checked_rnapower(depth = 20, n = 6, cv = 0.4, effect = 2, alpha = 0.05), 3), '\n')
cat('n for 80% power, 2-fold, CV=0.4:',
    ceiling(checked_rnapower(depth = 20, cv = 0.4, effect = 2, alpha = 0.05, power = 0.80)), 'per group\n')

# Effect / CV sensitivity: a SINGLE CV gives one curve, but real dispersion varies with mean.
for (cv in c(0.2, 0.4)) {
  n <- ceiling(checked_rnapower(depth = 20, cv = cv, effect = 1.5, alpha = 0.05, power = 0.80))
  cat(sprintf('CV %.1f -> n=%d per group for 1.5-fold at 80%% power\n', cv, n))
}

# ---------------------------------------------------------------------------
# 2. Simulation-based MARGINAL power at a target FDR (the honest default)
# ---------------------------------------------------------------------------
# Closed-form uses one CV; PROPER simulates from an empirical mean-dispersion trend and
# reports average power across the expression distribution; Actual FDR must pass below.
if (requireNamespace('PROPER', quietly = TRUE)) {
  library(PROPER)
  sim_opts <- RNAseq.SimOptions.2grp(ngenes = 20000, p.DE = 0.05,
                                     lOD = 'cheung', lBaselineExpr = 'cheung')
  sims <- runSims(Nreps = c(3, 5, 8, 12), sim.opts = sim_opts, nsims = 20, DEmethod = 'edgeR')
  powr <- comparePower(sims, alpha.type = 'fdr', alpha.nominal = 0.05,
                       stratify.by = 'expr', delta = log(1.5))   # delta is natural-log lfc in PROPER (not log2)
  fdr_gate <- assess_realized_fdr(summaryPower(powr), target_fdr = 0.05, tolerance = 0)
  print(fdr_gate)
  if (all(fdr_gate$accepted)) {
    cat('REALIZED_FDR_GATE=ACCEPTED\n')
  } else {
    cat('REALIZED_FDR_GATE=REJECTED; do not report/select these power values. Increase nsims, then refit the pilot dispersion/effect model or revise the DE method.\n')
  }
} else {
  cat('\nInstall PROPER for simulation-based marginal power (the reported figure).\n')
}

# ---------------------------------------------------------------------------
# 3. Depth vs replicates
# ---------------------------------------------------------------------------
# Past ~10-20M mapped reads, deeper sequencing adds little; replicates keep helping
# (Liu, Zhou & White 2014, Bioinformatics 30:301).
for (d in c(10, 20, 50, 100)) {
  cat(sprintf('Depth %3d, n=4: power = %.3f\n', d,
              checked_rnapower(depth = d, n = 4, cv = 0.4, effect = 2, alpha = 0.05)))
}
cat('# Adding a replicate beats doubling depth once depth is adequate.\n')
