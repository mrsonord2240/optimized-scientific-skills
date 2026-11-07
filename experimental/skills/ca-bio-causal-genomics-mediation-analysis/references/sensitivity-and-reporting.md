# Sensitivity and reporting
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
