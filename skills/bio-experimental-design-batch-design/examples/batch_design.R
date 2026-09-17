# Batch design: balanced assignment, hidden-batch detection, correction caveats
# Reference: sva 3.50+, designit 0.5+ | Verify API if version differs
#
# Demonstrates: the confounded-vs-balanced contrast, constrained sample-to-batch
# assignment, surrogate variable estimation for hidden batches, and the rule that a
# batch-"cleaned" matrix is for visualization only -- inference keeps batch in the model.

suppressPackageStartupMessages(library(sva))
set.seed(20260528)

# ---------------------------------------------------------------------------
# 1. Confounded vs balanced design
# ---------------------------------------------------------------------------
# BAD: condition aliased with batch -> non-identifiable, no correction can rescue it.
bad <- data.frame(condition = rep(c('ctrl', 'treat'), each = 12),
                  batch = rep(c('B1', 'B2'), each = 12))
cat('Confounded design (batch aliased with condition):\n'); print(table(bad$condition, bad$batch))

# GOOD: balance condition (and sex) across batches -> batch orthogonal, estimable.
samples <- data.frame(id = sprintf('S%02d', 1:24),
                     condition = rep(c('ctrl', 'treat'), each = 12),
                     sex = rep(c('M', 'F'), 12))
samples$batch <- NA_character_
for (cond in unique(samples$condition)) {
  idx <- which(samples$condition == cond)
  samples$batch[idx] <- sample(rep(paste0('B', 1:3), length.out = length(idx)))
}
cat('\nBalanced design (condition orthogonal to batch):\n'); print(table(samples$condition, samples$batch))

# Constrained assignment with designit (verify API vs installed vignette):
#   library(designit)
#   bc <- BatchContainer$new(dimensions = list(batch = 3, position = 8))
#   bc <- assign_in_order(bc, samples = samples)
#   bc <- optimize_design(bc, scoring = osat_score_generator(
#           batch_vars = 'batch', feature_vars = c('condition', 'sex')), max_iter = 10000)
#   assignment <- bc$get_samples()
#   # Verify before trusting the layout -- see SKILL.md "Verify the optimized layout"
#   tab <- table(assignment$condition, assignment$batch); print(tab)
#   stopifnot("condition confounded with batch" = all(tab > 0))

# Manual assignment above is a design-time balance choice too -- verify it the same way instead of
# only eyeballing the printed table:
tab_manual <- table(samples$condition, samples$batch)
stopifnot("condition is confounded with batch in the manual assignment" = all(tab_manual > 0))

# ---------------------------------------------------------------------------
# 2. Simulate a batch effect and detect hidden structure with SVA
# ---------------------------------------------------------------------------
n_genes <- 1000; n <- nrow(samples)
counts <- matrix(rnbinom(n_genes * n, mu = 100, size = 10), nrow = n_genes,
                 dimnames = list(paste0('Gene', 1:n_genes), samples$id))
batch_mult <- c(B1 = 1.0, B2 = 1.5, B3 = 0.7)              # multiplicative batch effect
counts <- sweep(counts, 2, batch_mult[samples$batch], '*')
counts[1:50, samples$condition == 'treat'] <- counts[1:50, samples$condition == 'treat'] * 2  # true DE
expr <- log2(counts + 1)

# Proteomics/metabolomics intensity matrices are rarely complete like this simulated RNA-seq count
# matrix -- inject missing values (MNAR-style: more likely at low abundance) to demonstrate that
# sva()/num.sv() require a complete matrix and to show the fix running on BOTH cases.
expr_na <- expr
set.seed(20260916)
na_prob <- pmin(1, pmax(0, (quantile(expr, 0.3) - expr) / quantile(expr, 0.3) * 0.6))
expr_na[matrix(runif(length(expr)) < na_prob, nrow(expr))] <- NA
cat(sprintf('\nInjected %d/%d NA cells (%.0f%%) into a copy of the matrix.\n',
            sum(is.na(expr_na)), length(expr_na), 100 * mean(is.na(expr_na))))

mod  <- model.matrix(~ condition, data = samples)
mod0 <- model.matrix(~ 1, data = samples)

run_sva_safely <- function(expr_normalized, mod, mod0, label) {
  n_bad <- sum(!is.finite(expr_normalized))
  if (n_bad > 0) {
    message(sprintf('[%s] %d/%d cells (%.0f%%) missing/non-finite; sva() requires a complete matrix.',
            label, n_bad, length(expr_normalized), 100 * n_bad / length(expr_normalized)))
    complete <- expr_normalized[stats::complete.cases(expr_normalized), , drop = FALSE]
    # Option A (used here): keep only complete-observation features -- biases toward abundant
    # features. Option B: impute first (normalization-qc: QRILC/missForest), then re-run.
    message(sprintf('[%s] Restricting to %d/%d complete-observation features (Option A).',
            label, nrow(complete), nrow(expr_normalized)))
    expr_normalized <- complete
  }
  stopifnot("expr_normalized still has non-finite values" = all(is.finite(expr_normalized)))
  n_sv <- num.sv(expr_normalized, mod, method = 'leek')
  cat(sprintf('[%s] rows used: %d | estimated hidden factors: %d\n', label, nrow(expr_normalized), n_sv))
  if (n_sv > 0) sva(expr_normalized, mod, mod0, n.sv = n_sv) else NULL
}

cat('\n--- SVA on the matrix as generated (no missing values) ---\n')
svobj_complete <- run_sva_safely(expr, mod, mod0, 'complete')

cat('\n--- SVA on the same matrix with injected missing values ---\n')
svobj_na <- run_sva_safely(expr_na, mod, mod0, 'NA-present')
# Correct use: add svobj$sv as covariates to the DE model (in differential-expression),
# NOT subtract them from `expr` before testing.

# ---------------------------------------------------------------------------
# 3. Cleaned matrix for VISUALIZATION ONLY (never feed into the hypothesis test)
# ---------------------------------------------------------------------------
# ComBat on log-normalized data; for integer counts use sva::ComBat_seq() instead.
expr_viz <- ComBat(dat = expr, batch = samples$batch, mod = mod)
pc_before <- cor(prcomp(t(expr))$x[, 1], as.numeric(factor(samples$batch)))
pc_after  <- cor(prcomp(t(expr_viz))$x[, 1], as.numeric(factor(samples$batch)))
cat(sprintf('PC1-vs-batch correlation: %.2f (raw) -> %.2f (ComBat, for plots only)\n',
            pc_before, pc_after))
# Inference rule (Nygaard 2016): keep batch in the model -- e.g. ~ condition + batch (+ SVs) --
# rather than testing on the ComBat-cleaned matrix, which understates residual variance.
