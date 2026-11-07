# Time-Varying Mediation and Double-ML (medDML)

### Time-Varying Mediation Methods

When the mediator is measured at multiple timepoints (or exposure varies over time), natural-effects estimands are not identified; switch to one of:

- g-formula (parametric or Monte Carlo): `gfoRmula::gformula_continuous_eof()`; `CMAverse::cmest(estimation='gformula')`
- g-estimation of a structural nested mean model: `gesttools::gestSingle()` / `gestMultiple()`
- Marginal structural model with stabilized IPTW: `ipw::ipwtm()` followed by `glm(..., weights=sw)`
- Sequential mediation for K timepoints: VanderWeele & Tchetgen Tchetgen 2017 JRSSB 79:917

Choose based on the experimental structure:
- >= 3 timepoints required for g-methods to identify time-varying indirect effects
- Longitudinal mediator measurement at EACH timepoint is required (not just baseline)
- MSM is preferred when treatment is binary and time-varying; g-formula when continuous
- Sequential mediation when the causal ordering of multiple mediators is known and stable across time

### Double-ML Doubly-Robust Mediation

**Goal:** Avoid model misspecification of both mediator and outcome models via cross-fitted ML nuisance estimators.

**Approach:** `causalweight::medDML` uses random forests (or other learners) with sample splitting to estimate nuisance parameters; final estimator is doubly robust.

```r
library(causalweight)

result_dml <- medDML(
  y=dat$outcome, d=dat$treatment, m=dat$mediator,
  x=as.matrix(dat[, covariates]),
  trim=0.05, order=1
)
# result_dml is a list of class "list" ($results, $ntrimmed) -- NOT a data.frame.
# $results is a 3x6 matrix: rows "effect"/"se"/"p-val" x columns
# total, dir.treat, dir.control, indir.treat, indir.control, Y(0,M(0))
# (verified on causalweight 1.1.4, 2026-09-17). Extract by name, not position:
total_effect <- result_dml$results['effect', 'total']
indirect_treat <- result_dml$results['effect', 'indir.treat']
indirect_control <- result_dml$results['effect', 'indir.control']
```

Reports direct (`dir.treat`/`dir.control`), indirect (`indir.treat`/`indir.control`, via mediator), and total effects with influence-function-based standard errors, plus `Y(0,M(0))` (baseline counterfactual mean). `dir.treat`/`indir.treat` and `dir.control`/`indir.control` are the effects evaluated with the mediator's treatment/control-arm distribution respectively (Farbmacher 2022's doubly-robust decomposition, analogous to `mediation::mediate()`'s treated/control ACME). Robust to non-linearity and interactions; assumes sequential ignorability still.
