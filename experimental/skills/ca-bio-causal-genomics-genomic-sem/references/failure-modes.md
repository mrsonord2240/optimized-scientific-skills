# Per-method failure modes (reference for genomic-sem SKILL.md)

Heywood cases, forcing a common factor on traits that do not share one, and a non-positive-definite V. Read it when `SKILL.md` "Reference Files" points here. The section names below refer to `SKILL.md`.

## Per-Method Failure Modes

### Heywood case (negative residual variance)

**Trigger:** A residual variance estimate is < 0, or a standardized loading exceeds 1.

**Mechanism:** Empirical underidentification; the genetic covariance matrix S is near-singular OR a trait has near-zero specific variance under the model. The maximum-likelihood / DWLS estimator runs past the boundary of the parameter space.

**Symptom:** `lavaan` warning "some estimated lv variances are negative" or "covariance matrix is not positive definite"; standardized loading > 1; non-convergence.

**Fix:** First, inspect the LDSC S matrix for genetic correlations near 1 (multicollinearity). Drop or merge near-identical traits. Second, constrain the offending residual variance to be non-negative in the lavaan syntax (`trait1 ~~ a*trait1; a > 0`). Third, verify the V matrix is positive definite via `chol(V_LD)`; if not, the bivariate LDSC inputs disagree on intercept signs and need re-munging. Never re-fit without diagnosing the cause.

### Trait inclusion under heterogeneous factor structure

**Trigger:** Forcing a common-factor model on traits that don't share a single latent factor.

**Mechanism:** When two or more traits load on a different factor than the rest, the single-factor model misfits. lavaan still returns parameter estimates but model fit is poor.

**Symptom:** CFI < 0.9; RMSEA > 0.08; some standardized loadings near 0 while others near 1; chi-square highly significant even after accounting for N.

**Fix:** Run ESEM first (`commonfactor` then `usermodel` with cross-loadings allowed) to discover structure. If two factors emerge, fit a two-factor `usermodel`. Drop traits with near-zero loadings on all factors. Document the model search.

### Non-positive-definite V_LD matrix

**Trigger:** LDSC inputs from different ancestry GWAS, or one trait with very low mean chi-square (< 1.02).

**Mechanism:** V is the sampling covariance of vech(S); when individual entries of S have huge SE relative to off-diagonal covariance, the resulting V is not positive definite (negative eigenvalues).

**Symptom:** `commonfactor()` errors with "matrix is not positive definite" before fitting; `eigen(LDSCoutput$V)$values` shows negative values.

**Fix:** Verify per-trait mean chi-square via `ldsc()` log; below 1.02, exclude that trait. Verify all GWAS are EUR ancestry (or match ancestry of LD scores). Apply nearest-PD smoothing via `Matrix::nearPD(V)$mat` ONLY as a last resort and document the approximation in methods.
