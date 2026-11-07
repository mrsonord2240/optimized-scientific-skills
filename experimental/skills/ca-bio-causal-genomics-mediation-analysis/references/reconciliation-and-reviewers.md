# Observational vs MR Reconciliation and Reviewer Pushback

## Reconciliation: Observational vs MR Mediation

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Observational ACME significant; MR-mediation null | Unmeasured M-Y confounding inflated observational estimate; OR weak IVs in MR | Re-run observational with `medsens()`; if rho_crit < 0.1, trust MR null |
| Observational ACME null; MR-mediation significant | Measurement error in M attenuated observational estimate | Trust MR (regression-dilution-free) IF instruments pass Steiger and pleiotropy tests (MR-Egger intercept, MR-PRESSO) |
| Both significant with same sign | Convergent evidence | High-confidence mediation; report effect size from the more-conservative estimate |
| Both significant with opposite signs | At least one is biased; revisit confounder structure and IV assumptions | Do not pool; investigate via cross-method sensitivity |

**Operational rule for high-stakes claims (clinical / drug-target mediation):** Require (1) significant observational ACME, (2) Imai rho_crit > 0.2 OR mediational E-value > 1.5, (3) directionally consistent MR-mediation result OR documented absence of valid instruments. Single-method mediation claims should be reported as exploratory.

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Sequential ignorability?" | Imai rho_crit reported via `medsens`; mediational E-value reported on the risk-ratio scale |
| "Exposure-induced confounder of M-Y?" | DAG drawn; if L present, switch to `CMAverse::cmest(estimation='msm')` for interventional indirect effect (Vansteelandt & Daniel 2017) |
| "Why this bootstrap method?" | BCa with sims=5000 for publication; percentile fallback when BCa fails to converge (acceleration estimate unstable at boundary) |
| "Why was MR-mediation not done?" | If valid IVs for E and M exist: two-step MR or MVMR-mediation done (see code below); if not, documented absence of trans-instruments |
| "Mediator measured with error?" | When mediator reliability r < 0.9, the indirect effect is attenuated. Regression calibration (Carroll 2006 Measurement Error in Nonlinear Models): replace the observed mediator with its conditional expectation given exposure and covariates. OR run a sensitivity analysis at fixed reliability r = 0.7 (Valeri L, Lin X & VanderWeele TJ 2014 Stat Med 33:4875) |
| "Why HIMA2 not BAMA?" | HIMA2 = frequentist + FDR control + faster; BAMA = Bayesian when prior information is available; sample-size justification given against simulation rule-of-thumb |
| "Proportion mediated unstable?" | When |total| < 2*SE(total), proportion-mediated CI is unreliable (denominator near zero); report indirect effect alone with absolute effect size |
