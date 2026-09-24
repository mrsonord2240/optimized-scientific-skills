# MR-Based Mediation (Two-Step and MVMR)

### Two-step MR instrument independence

**Trigger:** Same set of SNPs used as instruments for E in step 1 and for M in step 2.

**Mechanism:** If a SNP affects both E and M, the M-instrument violates exclusion restriction (the SNP-Y association is not exclusively through M). Estimates are biased toward the direct effect.

**Symptom:** Two-step MR shows large indirect effect; replacing M-instruments with non-overlapping SNPs makes it vanish.

**Fix:** Apply Steiger filter on the mediator: keep only SNPs where the SNP-M F-statistic exceeds SNP-E F-statistic (or where SNP explains more variance in M than E). For MVMR-mediation, require conditional F > 10 for each exposure independently (Sanderson 2019 IJE 48:713).

### MR-Mediation: Two-Step vs MVMR-Mediation

Decision tree:
- Independent instrument sets available for E and M -> two-step MR (Burgess 2015 IJE 44:484)
- E and M share instruments (common in cis-eQTL / cis-pQTL mediator cases) -> MVMR-mediation (Carter & Sanderson 2021 Eur J Epidemiol 36:465)
- Both feasible -> report both (triangulation)

Two-step code sketch:

```r
library(TwoSampleMR)
exp_E <- extract_instruments('ieu-a-2', clump=TRUE)
m_E <- extract_outcome_data(exp_E$SNP, 'ieu-b-30')
dat_EM <- harmonise_data(exp_E, m_E)
mr_EM <- mr(dat_EM)

exp_M <- extract_instruments('ieu-b-30', clump=TRUE)
exp_M_indep <- exp_M[!exp_M$SNP %in% exp_E$SNP, ]
exp_M_indep <- steiger_filtering(exp_M_indep)
out_M <- extract_outcome_data(exp_M_indep$SNP, 'ieu-a-7')
dat_MY <- harmonise_data(exp_M_indep, out_M)
mr_MY <- mr(dat_MY)
```

Indirect effect = beta_EM * beta_MY (product of coefficients). CI via delta method or parametric bootstrap of the joint (beta_EM, beta_MY) distribution. Steiger filter on M-instruments is mandatory to ensure the M -> Y direction (not Y -> M).

### MR-Mediation: Total Minus Direct via MVMR

**Goal:** Estimate the proportion of a genetic-instrument-identified causal effect that flows through a mediator, using independent IVs for E and (E + M).

**Approach:** Univariable MR for total E->Y; MVMR for direct E->Y conditional on M; indirect = total - direct via delta-method CI.

Runnable: `examples/mvmr_mediation.R` (`mr_ivw` total, `format_mvmr` + `strength_mvmr` + `ivw_mvmr` direct, indirect = total - direct with delta-method CI).

Require `fstat` conditional F > 10 for both E and M independently. If `fstat < 10`, use Q-statistic-adjusted IVW (`qhet_mvmr`) or report the result as weak-instrument-limited.
