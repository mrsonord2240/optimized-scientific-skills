# MVMR with conditional F

Read for multivariable MR (X1 adjusted for X2, mediation MR). The runnable code is `scripts/mvmr_conditional_f.R` (moved from SKILL.md "MVMR with Conditional F"); the prose below is verbatim.

## MVMR with Conditional F

**Goal:** Estimate the causal effect of exposure X1 on Y, adjusting for measured pleiotropy via X2.

**Approach:** Format exposures + outcome into MVMR object; compute conditional F per exposure (>10 required); run multivariable IVW; report Q_A heterogeneity.

```bash
Rscript scripts/mvmr_conditional_f.R --dat mvmr_input.tsv --exposures x1,x2 --gencov 0
```

`scripts/mvmr_conditional_f.R` (MVMR) reads a TSV with `SNP`, `beta.<x>`, `se.<x>` per exposure and `beta.y`, `se.y`, runs `format_mvmr()`, computes the per-exposure conditional F with `strength_mvmr()` (must be > 10 for EACH exposure, Sanderson 2019; total F is misleading), stops with an error at conditional F < 1 for any exposure (qhet_mvmr's own estimate is unreliable at this floor; can flip an exposure's sign; see the caveat below), then runs `ivw_mvmr()` and the Q_A heterogeneity test `pleiotropy_mvmr()`.

`gencov = 0` is valid ONLY if the exposure GWAS samples don't overlap; for overlapping exposures use the bivariate LDSC intercept matrix as `gencov`. If any conditional F < 10, the IVW point estimate is weak-IV-biased; the standard fallback is the Q-minimization estimator `qhet_mvmr(r_input, pcor, CI = TRUE, iterations = 1000)` (Sanderson 2021 Stat Med 40:5434), which minimizes Q-statistic heterogeneity rather than weighting by inverse variance.

**Caveat and guard (verified 2026-09-17, reproduced from an audit run; confidence-interval note added 2026-09-21):** qhet_mvmr's robustness to weak conditional instruments has a floor. At conditional F < 1 it does not just lose precision -- it can flip the sign of an exposure's estimate. On synthetic instruments with planted direct effects 0.30 / -0.10 and conditional F = 0.87 / 0.78, MVMR-IVW stayed close to truth (0.298 / -0.099) but qhet_mvmr gave 0.236 / **+0.049** -- the second exposure's sign flipped. Guard: below conditional F = 1, do not trust qhet_mvmr's point estimate as a correction for exposure 1 or 2; report the weak-IV-biased MVMR-IVW estimate with that caveat instead, or acquire stronger/less-correlated instruments before drawing a directional conclusion. The code above raises an error at that floor rather than silently returning a fallback estimate that may be flipped. Above the floor the estimate is correctly signed but imprecise until conditional F nears 10 (a re-audit run at conditional F ~2.2, planted 0.30 / -0.10: qhet_mvmr gave 0.44 [0.12, 0.88] / -0.26 [-0.67, 0.16]), so report it with its confidence interval, never as a precise correction.
