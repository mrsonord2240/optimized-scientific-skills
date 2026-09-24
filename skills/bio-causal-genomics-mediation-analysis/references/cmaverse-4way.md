# CMAverse 4-Way Decomposition and Exposure-Induced Confounders

## 4-Way Decomposition Framework

VanderWeele 4-way explicitly separates effects from exposure-mediator interaction:

```
Total Effect = CDE + INTref + INTmed + PIE
```

| Component | Meaning | Active when |
|-----------|---------|-------------|
| CDE | Controlled direct effect (with mediator fixed at reference level) | E directly affects Y |
| INTref | Interaction-reference -- needs interaction AND exposure | E*M interaction with mediator at reference |
| INTmed | Mediated interaction -- needs interaction AND exposure AND mediation | E shifts M which then interacts with E |
| PIE | Pure indirect effect (older "mediation" quantity) | E shifts M which shifts Y additively |

Without an exposure-mediator interaction term, INTref = INTmed = 0 and the decomposition collapses to CDE + PIE (= ADE + ACME). With interaction present, traditional ACME mixes PIE and INTmed; the 4-way separation is the only framework that disentangles them. Most epidemiology applications include the interaction term and report all four components (Valeri & VanderWeele 2013 Psychol Methods).

### Exposure-induced M-Y confounder

**Trigger:** A covariate L sits between E and Y, AND is affected by E, AND confounds M-Y.

**Mechanism:** Standard regression-based mediation cannot adjust for L without blocking part of the indirect effect (collider stratification bias). Adjusting biases CDE; not adjusting biases ACME.

**Symptom:** Sensitivity to confounder set; ACME flips sign when L is added vs removed.

**Operational identification:** From the DAG, L is a covariate of M and Y that is also affected by E. VanderWeele TJ, Vansteelandt S & Robins JM 2014 (Epidemiology 25:300) give the criterion: if L is adjusted, part of the indirect E -> L -> M -> Y pathway is blocked; if L is not adjusted, L confounds the M-Y leg. Both are wrong under natural-effects; the natural indirect effect is simply not identified.

**Fix:** Switch to **interventional indirect effects** (Vansteelandt & Daniel 2017 Epidemiology 28:258), NOT natural indirect effects. Use `CMAverse::cmest(estimation='msm')` with stabilized inverse-probability weights (yields the randomized-interventional analogue), `gfoRmula` (parametric g-formula), or randomized/interventional indirect effects (Lin SH & VanderWeele TJ 2017 J Causal Inference 5:20150027). The interventional indirect is identified under weaker assumptions than the natural indirect.

### 4-Way Decomposition with Exposure-Mediator Interaction

**Goal:** Separate CDE, PIE, INTref, INTmed (continuous outcome) -- or their excess-relative-risk equivalents ERcde/ERpnie/ERintref/ERintmed (binary/survival outcome, shown below) -- when exposure-mediator interaction is biologically plausible (e.g., gene-environment interaction modifying mediator effect).

**Approach:** Use CMAverse regression-based estimator with `EMint=TRUE`; bootstrap CIs.

Runnable: `examples/cmaverse_4way.R` (`cmest(model='rb', EMint=TRUE, mreg=list('linear'), yreg='logistic', astar=0, a=1, mval=list(0), estimation='paramfunc', inference='bootstrap', nboot=1000)`, plus the mediational E-value and probit `medsens()`).

CMAverse reports the 4-way decomposition (Vanderweele 2014): for continuous outcomes the components are `cde`, `intref`, `intmed`, `pnie` (or `pie`), `te`, `pm`; for non-continuous outcomes (logistic / Cox / Poisson) the ratio effects `Rcde`, `Rpnde`, `Rtnde`, `Rpnie`, `Rtnie`, `Rte` are reported ALONGSIDE an excess-relative-risk decomposition with an `ER` prefix -- `ERcde`, `ERintref`, `ERintmed`, `ERpnie` (plus a `(prop)` share for each) -- **not** the bare `intref`/`intmed` names, which belong only to the continuous-outcome case. When `EMint=TRUE`, `pm`, `int`, `pe` are also included. Verified column set on a logistic-outcome, `EMint=TRUE` fit (CMAverse 0.1.0, 2026-09-17): `Rcde Rpnde Rtnde Rpnie Rtnie Rte ERcde ERintref ERintmed ERpnie ERcde(prop) ERintref(prop) ERintmed(prop) ERpnie(prop) pm int pe`. Verify column names in the installed CMAverse version with `summary(result)$summarydf` (not `$results`, which does not exist on the summary object), since naming has evolved.
