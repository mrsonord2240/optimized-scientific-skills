---
name: bio-causal-genomics-mediation-analysis
category: Data Analysis
description: Decompose total effects into direct and indirect paths through mediators using mediation, CMAverse 4-way, HIMA/HIMA2 high-dimensional, BAMA, two-step / MVMR mediation, or double-ML medDML. Use when testing whether a molecular phenotype (expression, methylation, protein) mediates a treatment-outcome relationship, decomposing exposure-mediator interaction via VanderWeele 4-way, screening high-dimensional EWAS mediators, or running MR-based mediation when sequential ignorability is implausible.
tool_type: mixed
primary_tool: mediation
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: R 4.3+, mediation 4.5.0+, CMAverse 0.1.0+ (GitHub `BS1125/CMAverse`), HIMA >= 2.3.0 (CRAN; checked on HIMA 2.3.4), bama 1.3+, causalweight 1.0.5+ (medDML), MVMR 0.4+, TwoSampleMR 0.6+, EValue 4.1+, gesttools 1.3+, ipw 1.0.11+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- HIMA must be pinned at `>= 2.3.0` for the code patterns below; the formula interface `hima(formula, data.pheno, data.M, mediator.type, penalty, ...)` was introduced in 2.3.0. On HIMA 2.2.x the classic engine is the top-level `hima(X, Y, M, COV.XM=, COV.MY=, Y.family=, penalty=)` (positional order X, Y, M) with no formula interface; HIMA 2.2.x is NOT API-compatible with the examples here.
- In 2.3+ the classic engine is `hima_classic(X, M, Y, COV.XM=, COV.MY=, Y.type=)` (positional order X, M, Y; outcome type via `Y.type`, not `Y.family`); use it only to reproduce the original 2016-2021 SIS+penalty pipelines.

If code throws an error, introspect the installed package (`?hima`, `args(cmest)`) and adapt the example to match the actual API rather than retrying.

# Mediation Analysis

**"Does expression of GENE_X mediate the SNP-to-disease effect?"** -> Decompose the total effect of a treatment (genotype, exposure) on an outcome into direct and indirect paths through one or more mediators, with explicit handling of exposure-mediator interaction, sensitivity to unmeasured confounding, and high-dimensional mediator screening.

- R (single-mediator, observational, sequential-ignorability assumed): `mediation::mediate(med_model, out_model, treat='X', mediator='M', boot=TRUE, sims=5000)`
- R (4-way decomposition with exposure-mediator interaction): `CMAverse::cmest(...EMint=TRUE, estimation='paramfunc', inference='bootstrap', nboot=1000)`
- R (high-dimensional / EWAS mediators): `HIMA::hima(Y ~ X + covariates, data.pheno, data.M, mediator.type='gaussian', penalty='DBlasso')` (modern v2.3+ formula interface)
- R (MR-based mediation): two-step `TwoSampleMR` with independent instruments OR `MVMR::ivw_mvmr` for joint direct effect
- R (doubly-robust double-ML): `causalweight::medDML(y, d, m, x)`

Sequential ignorability (no unmeasured confounder of treatment-mediator, mediator-outcome, treatment-outcome) is the single load-bearing assumption of observational mediation and is fundamentally untestable. Every report should include a sensitivity result (Imai's rho via `medsens()` or a mediational E-value).

## Scope

ACME/CDE/PIE/indirect-effect estimates are population-level causal-inference quantities from a fitted model or GWAS under an assumed DAG. They do not license an individual patient's treatment decision -- decline requests to turn a population-level mediation result into an individual dosing or treatment recommendation, and redirect to the treating clinician, who has clinical context this Skill does not.

Individual-level genotype, expression and methylation data used as mediators are PHI-adjacent. De-identify before analysis, and avoid echoing raw per-subject values in error messages, logs or intermediate output.

## Algorithmic Taxonomy

| Method | Framework | Handles E-M interaction | High-D mediators | Min n | Fails when |
|--------|-----------|--------------------------|------------------|-------|------------|
| Baron-Kenny (1986) | Additive regression-based product/difference | No | No | ~100 | Any non-linearity, interaction, or binary outcome; deprecated for causal inference |
| Imai mediation R (Imai 2010 Psychol Methods) | Counterfactual ACME/ADE with bootstrap | Yes (via interaction term in outcome model) | No | ~200 | Sequential ignorability violated; exposure-induced M-Y confounder; rare binary outcome with logistic outcome model |
| VanderWeele 4-way (VanderWeele 2014 Epidemiology 25:749; 2015 OUP) | CDE + PIE + INTref + INTmed decomposition | Native | No | ~300 | Without interaction term reduces to standard mediation; binary outcome needs rare-disease assumption |
| CMAverse (Shi 2021 Epidemiology 32:e20) | 6 estimators: regression (`rb`), weighting (`wb`), IORW (`iorw`), natural effect models (`ne`), MSM (`msm`), g-formula (`gformula`) | Yes | No (single M, or M-vector) | ~300 | Estimator-specific; `wb` fails with rare exposure; `msm` needs censoring weights for survival |
| HIMA1 / hima_classic (Zhang 2016 Bioinformatics 32:3150) | SIS screen by beta (M->Y) + MCP penalty | No | Yes (up to ~10k) | ~150 + p>>n | Misses mediators with strong alpha and weak beta; screening-step false-negatives |
| HIMA2 / hima (Perera 2022 BMC Bioinformatics 23:296) | SIS screen by alpha*beta (indirect effect) + de-biased Lasso (DBlasso) | No (linear by default) | Yes | ~150 | Outcome family limited (gaussian/binomial); HIMA-Cox for survival; HIMA-Pois for count |
| HILAMA (Wang et al 2025) | High-D mediation with latent confounders | No | Yes (>= 100k) | ~500 | Newer; benchmarks evolving; requires latent-factor specification |
| BAMA (Song 2020 Biometrics) | Bayesian high-D continuous shrinkage | No | Yes (~5k) | ~200 | Slow MCMC; prior sensitivity for very weak mediators |
| Two-step MR / network MR (Burgess 2015 IJE 44:484) | IV-based at each step with INDEPENDENT instruments | Implicit (no interaction modeling) | One mediator at a time | Large summary-stat samples | Same SNP used for E and M (violates exclusion); horizontal pleiotropy; Steiger reversal of M->E direction |
| MVMR-mediation (Carter & Sanderson 2021 Eur J Epidemiol 36:465-478) | Total minus direct via MVMR | Implicit | Single mediator | Large GWAS samples for both E and M | Conditional F < 10 for either exposure; correlated instruments |
| medDML (Farbmacher 2022 Econometrics J 25:277) | Double-debiased ML, doubly-robust | Limited (depends on learner) | Moderate (sparsity-friendly) | ~500 | Severe overlap violations; cross-fitting variance with small n |

Methodology evolves; verify against the current CMAverse vignette and the Steen / Vansteelandt natural-effects-model literature before locking analytic choices. Difference-in-coefficients and product-of-coefficients give identical estimates in fully linear-Gaussian models but DIVERGE for any non-linear outcome model (logistic, Cox, Poisson); the counterfactual ACME from `mediation::mediate()` is the correct quantity for non-linear outcomes.

## Decision Tree by Scenario

| Scenario | Recommended pipeline |
|----------|---------------------|
| Observational, single measured mediator, no plausible E-M interaction, continuous outcome | `mediation::mediate()` with `boot=TRUE, sims=5000`; always run `medsens()` |
| Observational, single mediator, suspected E-M interaction, any outcome family | `CMAverse::cmest(..., EMint=TRUE)` -> read `cde`/`intref`/`intmed`/`pnie` (continuous outcome) or `ERcde`/`ERintref`/`ERintmed`/`ERpnie` (ratio-scale, e.g. logistic/Cox); verify with `summary(result)$summarydf`; see `references/cmaverse-4way.md` |
| Observational, BINARY outcome, rare disease (< 10%) | `cmest(yreg='logistic', EMint=TRUE, casecontrol=FALSE)` -- OR-based 4-way decomposition is valid under rare-disease |
| Observational, survival outcome | `cmest(yreg='coxph')` OR `HIMA::hima_cox` for high-D; report HRs |
| High-D mediators (EWAS, transcriptome-wide), continuous outcome | `HIMA::hima(formula, data.pheno, data.M, mediator.type='gaussian', penalty='DBlasso')`; report `sigcut` (FDR threshold, default 0.05); see `references/hima-ewas.md` |
| High-D mediators with latent confounding (very-high-D EWAS) | `HILAMA` (2025) |
| High-D mediators with Bayesian shrinkage (small n, ~5k features) | `bama::bama()` |
| Strong genetic IVs for exposure available, single mediator with own IVs | Two-step MR with independent instruments + Steiger filter on mediator; see `references/mr-mediation.md` |
| Both E and M have IVs but instruments are weak / correlated | MVMR-mediation with conditional F > 10 each; see `references/mr-mediation.md` |
| Observational with rich confounder set, want doubly-robust estimate | `causalweight::medDML` (double-debiased ML); see `references/time-varying-and-dml.md` |
| Longitudinal with time-varying confounding | g-formula via `CMAverse::cmest(estimation='gformula')` OR `gfoRmula` package; see `references/time-varying-and-dml.md` |
| Exposure-induced confounder of M-Y exists | Interventional indirect effects (Vansteelandt & Daniel 2017); `CMAverse::cmest(estimation='msm')`; see `references/cmaverse-4way.md` |

BAMA vs HIMA2 for high-D mediators: prefer BAMA over HIMA2 when (1) strong prior information on mediator effects is available, (2) the candidate panel is moderately sized (~5k mediators), and (3) the compute budget allows 1-6h MCMC; otherwise HIMA2 is faster with comparable FDR control.

## Sequential Ignorability and Why It Always Needs Sensitivity

Observational mediation requires three no-unmeasured-confounding assumptions. The third (M-Y unmeasured confounder, after conditioning on E) is the most common violator in genomic mediation because biological confounders (cell composition, batch effects, technical mediators) frequently affect both M and Y.

### Sequential ignorability untestable

**Trigger:** Always, by design.

**Mechanism:** No statistical test can detect an unmeasured confounder of M-Y. Bootstrap CIs assume the assumption holds; they do NOT propagate uncertainty about it.

**Symptom:** Significant ACME with no sensitivity reported -> reviewer rejects.

**Fix:** Report at least one of:
- Imai's rho sensitivity: `medsens(med_result, rho.by=0.05, sims=1000)`; the critical rho where ACME crosses 0; |rho_crit| > 0.3 is "reasonably robust" (Imai 2010), |rho_crit| < 0.1 is highly sensitive.
- Mediational E-value (Smith & VanderWeele 2019 Epidemiology 30:835): minimum risk-ratio strength of an unmeasured confounder needed to nullify the observed indirect effect; computed via `EValue::evalues.OLS()` for linear outcomes or by-hand from ACME risk ratio bounds.
- Reporting BOTH rho-based and E-value sensitivity is standard for high-stakes claims.

### Methods-Section Defense of Sequential Ignorability

Template sentence for the methods write-up: "We assumed sequential ignorability conditional on {age, sex, ancestry PCs, cell composition, batch, smoking}. Robustness was assessed via Imai rho_crit at the ACME contrast (`medsens`, sims = 1000) and the mediational E-value on the risk-ratio scale (Smith & VanderWeele 2019 Epidemiology 30:835)."

Quantitative interpretation thresholds:
- rho_crit: |rho_crit| > 0.3 robust; 0.1-0.3 moderately sensitive; < 0.1 highly sensitive (Imai 2010).
- E-value: E > 2 robust to plausible biological confounding; 1.5-2 moderate; < 1.5 fragile (Smith & VanderWeele 2019).

For high-stakes claims (clinical, drug-target, regulatory submissions) report BOTH rho_crit and the mediational E-value; for exploratory work either alone suffices.

### Bootstrap iterations too low

**Trigger:** `sims=100` or `sims=500` in early exploration left in for the final report.

**Mechanism:** ACME CIs from bootstrap have Monte-Carlo error that scales as 1/sqrt(sims); at sims=500 the 95% CI bounds have ~5% MC noise, enough to flip the conclusion at the boundary.

**Symptom:** Re-running `mediate()` with a different `set.seed()` gives substantially different CI bounds.

**Fix:** `sims=1000` minimum for any reported result; `sims=5000` for publication; `sims=10000` if proximity to zero matters. BCa CIs (`boot.ci.type='bca'` in `mediate()` -- lowercase `'bca'`, not `'BCa'`) are slightly more accurate than percentile CIs near zero but require more sims for stability.

### Difference vs product of coefficients diverge for non-linear outcomes

**Trigger:** Binary or survival outcome modeled with logistic / Cox.

**Mechanism:** Difference = total - direct; product = alpha * beta. Equivalent under linear-Gaussian; diverge under any link function. Counterfactual ACME from `mediation::mediate()` is the correct quantity; hand-computed product-of-coefficients on logistic output is biased except under rare-disease.

**Fix:** Report only counterfactual ACME (Imai or CMAverse). For OR-based 4-way decomposition on rare outcomes (<= 10%), Valeri & VanderWeele 2013 formulas apply; for common outcomes use risk-ratio scale or marginal effects rather than ORs.

### Required Reporting for Publication

| Component | Required |
|-----------|----------|
| ACME estimate + 95% CI (BCa preferred) | Yes |
| ADE + 95% CI | Yes |
| Total effect | Yes |
| Proportion mediated | Yes when total > effect-size threshold |
| Bootstrap method + sims | percentile / BCa; min 1000, recommend 5000 |
| Sequential ignorability sensitivity | rho_crit (medsens) OR mediational E-value |
| Exposure-mediator interaction test | Coefficient + p; 4-way decomposition if significant |
| Confounder set justification | DAG description |
| Mediator measurement reliability | Cite |
| Sample size + missing-data handling | Yes |
| Mediator / exposure scale | Standardized? log? raw? |

Reference: AGReMA guideline (Lee H et al 2021 JAMA 326:1045) and MacKinnon 2008 Introduction to Statistical Mediation Analysis.

## Quantitative Thresholds

| Quantity | Threshold | Source / Rationale |
|----------|-----------|--------------------|
| `sims` (bootstrap iterations) | >= 1000 exploratory, >= 5000 publication | Imai 2010; MC error scales 1/sqrt(sims) |
| Proportion mediated -- meaningful | > 0.2 | Convention; weak guideline only -- effect size in absolute terms matters more (MacKinnon 2008) |
| Proportion mediated -- "most of the effect" | > 0.5-0.8 | Convention |
| Imai rho_crit -- robust | > 0.3 | Imai 2010 Psychol Methods 15:309 |
| Imai rho_crit -- sensitive | < 0.1 | Same |
| Mediational E-value -- robust | > 2.0 (working convention; the original Smith & VanderWeele 2019 E-value framework does not prescribe a specific cutoff -- magnitude is context-dependent) | Smith & VanderWeele 2019 Epidemiology 30:835 |
| HIMA FDR cutoff | BH FDR < 0.05 | Default; report q-values not raw p |
| MVMR conditional F per exposure | > 10 each | Sanderson 2019 IJE 48:713 |
| Two-step MR -- F for both stages | > 10 each | Burgess weak-instrument convention |
| Sample size -- single-mediator (Imai) | >= 200 for stable bootstrap | Simulation rule-of-thumb |
| Sample size -- HIMA EWAS | >= 150 with p_mediators up to ~10k | Zhang 2016 simulations |
| Rare-outcome cutoff for OR-based 4-way | outcome prevalence <= 10% | Valeri & VanderWeele 2013 |

## Working Code Patterns

### Single-Mediator Observational with Sensitivity

**Goal:** Decompose a genotype-disease effect via measured molecular mediator with explicit sensitivity to unmeasured M-Y confounding.

**Approach:** Fit mediator and outcome models, bootstrap ACME/ADE, run `medsens()` for Imai rho-based sensitivity.

Runnable: `examples/eqtl_mediation.R` (mediator + outcome models, `mediate(..., boot=TRUE, sims=1000)`, multi-gene loop; use `sims=5000` and `boot.ci.type='bca'` for publication) and `examples/sensitivity_analysis.R` (`medsens(..., rho.by=0.05, effect.type='indirect')`, which needs a **probit** outcome model -- see Common Errors).

`d0`, `z0`, `n0`, `tau.coef` slots return ACME, ADE, proportion mediated, total effect.

### Mediational E-Value for Sensitivity

**Goal:** Report the minimum strength of an unmeasured M-Y confounder required to nullify the observed indirect effect.

**Approach:** Convert ACME and its CI to a risk-ratio scale, then apply VanderWeele E-value formula.

```r
library(EValue)

acme_rr <- exp(med_result$d0)
acme_lower_rr <- exp(med_result$d0.ci[1])
evalues.RR(acme_rr, lo=acme_lower_rr)   # omit `hi`: its default is NA; `hi=NULL` crashes (`if (est < true & !is.na(hi))`, EValue 4.1.4)
```

For binary outcomes, convert ACME on probability scale to RR; for continuous, use `evalues.OLS()` with the standardized indirect effect. E-value > 2 indicates a confounder would need >2-fold associations with both M and Y to nullify the indirect effect (Smith & VanderWeele 2019). By hand: convert ACME to a risk-ratio bound (`acme_rr = exp(ACME)` on the log scale for continuous outcomes, or VanderWeele's marginal RR conversion for binary), then `E = RR + sqrt(RR * (RR - 1))`; apply the same formula to the CI bound closer to the null for the E-value of the CI. `EValue::evalues.OLS()` automates this for linear outcomes.

## Reference Files

Read the file for the method in use; everything a request always needs (scope, taxonomy, decision tree, sensitivity, thresholds, install, Common Errors) stays in this file.

| File | Read when |
|------|-----------|
| `references/cmaverse-4way.md` | 4-way CDE/INTref/INTmed/PIE decomposition, `cmest()` output names, or an exposure-induced M-Y confounder (interventional indirect effects) |
| `references/hima-ewas.md` | HIMA/HIMA2 on high-dimensional mediators: `data.pheno` covariate errors, `mediator.type`, the `hima()` code pattern, correlated mediators |
| `references/mr-mediation.md` | Two-step MR (instrument independence, Steiger) or total-minus-direct MVMR mediation |
| `references/time-varying-and-dml.md` | Longitudinal mediator/exposure (g-formula, SNMM, MSM) or `causalweight::medDML` double-ML |
| `references/reconciliation-and-reviewers.md` | Observational ACME and MR-mediation disagree, or drafting responses to reviewers |

Runnable code: `scripts/hima_ewas.R` (HIMA pipeline); `examples/` (`eqtl_mediation.R`, `sensitivity_analysis.R`, `cmaverse_4way.R`, `mvmr_mediation.R`).

## Tool Install Notes

| Package | Source | Notes | Compute time |
|---------|--------|-------|--------------|
| mediation | CRAN | `install.packages('mediation')`; actively maintained (Imai group) | Single mediator, 5000 bootstrap sims: minutes on a laptop |
| CMAverse | GitHub | `remotes::install_github('BS1125/CMAverse')`; NOT on CRAN; 6 estimators in one interface | Minutes for `nboot=1000` on a single dataset |
| HIMA | CRAN | `BiocManager::install('qvalue')` then `install.packages('HIMA')`; back on CRAN with a lighter dependency list (`ncvreg`, `glmnet`, Bioconductor `qvalue`), no `scalreg` needed; confirmed on HIMA 2.3.4 (checked 2026-09-17) | ~500k CpGs, n=500: 30-60 min with `parallel=TRUE, ncore=8` |
| bama | CRAN | `install.packages('bama')`; Bayesian; slow MCMC | 1-6 hours depending on chain length |
| causalweight | CRAN | `install.packages('causalweight')`; medDML for double-ML mediation | Cross-fitted random forests: 10-30 min for n=2000 |
| EValue | CRAN | `install.packages('EValue')`; for mediational E-values | Seconds |
| TwoSampleMR | r-universe | See causal-genomics/mendelian-randomization for setup | Depends on instrument count and OpenGWAS API latency |
| MVMR | r-universe | `remotes::install_github('WSpiller/MVMR')`; for MVMR-mediation | Seconds to minutes |
| gfoRmula | CRAN | For longitudinal / time-varying confounders | Minutes to tens of minutes, scales with timepoints |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `Error in storage.mode(x) <- "double"` inside `hima()` | NA in `data.pheno` columns referenced by formula, or unconverted factors | `na.omit(data.pheno)` first; ensure all RHS vars in formula are numeric or dummy-coded via `model.matrix()` (HIMA 2.3.4 rejects `factor` columns; see `references/hima-ewas.md`) |
| ACME significant, ADE significant, total NOT significant | Suppression / inconsistent mediation | Report transparently; effect partitioning can exceed total in suppression |
| `medsens()` errors on glm outcome | `medsens` requires linear OR probit (not logit) outcome | Refit outcome as `glm(..., family=binomial(link='probit'))` |
| `mediate()` runs forever with binary outcome | `sims=5000` with bootstrap and small n | Use `sims=1000` exploratory; verify model converges first; consider parallel via `parallel='multicore'` |
| CMAverse `cmest()` reports NaN for `pm` | Total effect crosses zero -> proportion ill-defined | Report ACME and TE separately; pm is unstable when |TE| is small |
| Different ACME between `mediation` and CMAverse `rb` | Default `astar/a` levels differ; binary mediator handled differently | Set `astar=0, a=1` explicitly; for binary mediator pass `mval=list(0)` |
| HIMA returns zero significant mediators | Screening too aggressive; or no true mediators | Try `topN=2*sqrt(n)` instead of default; verify with permutation null |
| `nrow(result)`/`rownames(result)` is `NULL` after `hima()` (no error thrown) | `hima()` returns a list of class `"hima"`, not a data.frame; these accessors fail silently | Use `result$ID` (mediator names), `length(result$ID)` (count); index `result$alpha`, `result$beta`, `result$rimp` the same way |
| Two-step MR shows indirect > total | Steiger reversal: M actually causes E; or pleiotropic SNPs | Run MR-Steiger filter; use MR-PRESSO for pleiotropy |
| `medDML` trim removes most data | Severe positivity violation -- few units with overlapping treatment/mediator distributions | Tighten covariate set; check propensity score distributions |
| `medDML` fails with `Error ... subscript out of bounds` | Internal cross-fitted Lasso step (`hdm::rlassologit`) needs a covariate matrix with named columns and enough independent variation; too-few observations or near-collinear/unnamed `x` columns break it (verified: n=150, 2 unnamed near-collinear columns crashes; n=800 with named `age`/`sex`/`bmi` columns converges cleanly, causalweight 1.1.4, 2026-09-17) | Use `x=as.matrix(dat[, covariates])` with named columns (not a bare `cbind()` of unnamed vectors), n >= ~500, and drop near-collinear covariates |

## References

- Baron RM, Kenny DA 1986 J Pers Soc Psychol 51:1173 (original product-of-coefficients)
- Imai K, Keele L, Tingley D 2010 Psychol Methods 15:309 (counterfactual mediation, sequential ignorability)
- VanderWeele TJ 2014 Epidemiology 25:749 (4-way decomposition)
- VanderWeele TJ 2015 Explanation in Causal Inference (OUP) -- canonical textbook
- Valeri L, VanderWeele TJ 2013 Psychol Methods 18:137 (binary outcomes; rare-disease 4-way)
- Vansteelandt S, Daniel RM 2017 Epidemiology 28:258 (interventional / randomized indirect effects)
- Shi B et al 2021 Epidemiology 32:e20 (CMAverse package; 6 estimators)
- Zhang H et al 2016 Bioinformatics 32:3150 (HIMA original)
- Perera C et al 2022 BMC Bioinformatics 23:296 (HIMA2 alpha-beta screening)
- Song Y et al 2020 Biometrics 76:700 (BAMA Bayesian high-D mediation)
- Burgess S et al 2015 IJE 44:484 (network / two-step MR)
- Sanderson E et al 2019 IJE 48:713 (MVMR conditional F-statistic)
- Carter AR & Sanderson E 2021 Eur J Epidemiol 36:465 (MVMR-mediation)
- Farbmacher H et al 2022 Econometrics J 25:277 (medDML / double-ML mediation)
- Smith LH, VanderWeele TJ 2019 Epidemiology 30:835 (mediational E-value)

## Related Skills

- causal-genomics/mendelian-randomization - IV-based causal inference; foundation for MR-mediation
- causal-genomics/pleiotropy-detection - MR-mediation instrument validity hinges on pleiotropy diagnostics (MR-Egger intercept, MR-PRESSO)
- causal-genomics/colocalization-analysis - Confirm shared causal variant before causal mediation
- causal-genomics/fine-mapping - Identify the causal variant driving the exposure
- methylation-analysis/differential-cpg-testing - Per-CpG inputs for HIMA EWAS mediation
- differential-expression/deseq2-basics - Expression inputs for eQTL mediation
- multi-omics-integration/mofa-integration - Multi-layer mediator construction
- population-genetics/association-testing - GWAS summary statistics for MR-mediation
- clinical-biostatistics/effect-measures - Risk-ratio / odds-ratio scales for binary outcomes
- machine-learning/model-validation - Cross-fitting and sample splitting for medDML
