# Response-Adaptive Randomisation, Platforms, and Dose Finding

Use this reference for multi-arm allocation changes, platform entry/exit,
treatment graduation, and Phase 1 dose-finding decisions.

## Response-adaptive randomisation

RAR changes allocation probabilities as outcomes accumulate. Its patient-
welfare argument is strongest in multi-arm (at least three-arm), rare-disease,
or biomarker-stratified settings and weakest in two-arm confirmatory trials.
Time trends can be confounded with allocation changes, and patients can mistake
adaptive allocation for individualised treatment.

For any RAR design, pre-specify:

- the outcome model and update cadence;
- burn-in and minimum allocation probabilities;
- stratification and calendar-time covariates;
- delayed-outcome handling and data cutoffs;
- randomisation and primary-analysis weights;
- stopping, graduation, and arm-entry rules;
- simulations with secular drift, site learning, recruitment shifts, missing
  outcomes, and model misspecification.

The modern consensus reflected by Robertson et al. (2023) does not justify RAR
merely because software can perform it. For a two-arm confirmatory trial, prefer
group-sequential efficacy and futility boundaries unless a compelling,
simulated case establishes otherwise.

## Bayesian platform trials

I-SPY 2 is a multi-arm, biomarker-stratified Bayesian platform using RAR and a
graduation rule based on posterior predictive probability of success in a
future trial. Its commonly described criterion is at least 85% predicted
probability of success in a 300-patient Phase 3 trial. Treat those values as an
I-SPY example, not general defaults.

GBM AGILE illustrates a registrational global oncology platform. REMAP-CAP
illustrates a Bayesian factorial multi-domain design. Each requires explicit
domain/arm entry and exit rules, shared-control handling, multiplicity and
time-trend evaluation, and versioned simulation.

Adaptive futility closure is statistically simpler than graduating the most
promising arm. Graduation selects on a favourable interim estimate, so report a
bias-adjusted estimator, such as conditional maximum likelihood, alongside the
raw estimate.

## Phase I dose finding

| Design | Core idea | Main advantage | Main risk |
|---|---|---|---|
| CRM | Parametric dose-toxicity model updated by cohort | Efficient when the skeleton is calibrated | Skeleton misspecification can redirect the trial |
| EWOC | CRM-like model with an overdose-control constraint | Explicit overdose protection | Still depends on model calibration |
| mTPI | Beta-binomial unit-probability-mass rule | Pre-tabulated decisions | Original formulation has an interval-width/Ockham bias |
| mTPI-2 / Keyboard | Equal-width interval refinement | Removes the original mTPI interval bias | Requires calibrated operating characteristics |
| BOIN | Optimised interval escalation/de-escalation boundaries | Transparent bedside table; near-CRM performance | Target, cohort, and stopping rules still require simulation |

BOIN was accepted into the FDA Fit-for-Purpose program in December 2021. That
does not make every BOIN design acceptable. Generate and lock the escalation
table, define overdose and early-stop rules, and simulate the true dose-
toxicity curves most relevant to the program.

For CRM, calibrate the skeleton, for example with the Lee-Cheung indifference-
interval method. Compare MTD selection, overdose allocation, early stopping,
and expected sample size against BOIN or another justified comparator.

Project Optimus displaced the simple MTD-and-go pattern in oncology. After dose
escalation, use randomised multi-dose comparison where appropriate, integrating
exposure, safety, tolerability, pharmacodynamics, and preliminary activity.

## Basket and borrowing designs

EXNEX or robust meta-analytic predictive priors let strata borrow while
retaining a non-exchangeable component that permits an outlying stratum to
detach. Pre-specify exchangeability weights and test sensitivity to prior-data
conflict. Pediatric power-prior discount values such as 0.3-0.6 are examples to
stress-test, not universal defaults.

The packaged R script demonstrates only a single historical-data MAP prior. It
does not implement stratum-level EX/NEX components or conflict-driven
detachment and must not be presented as an EXNEX basket analysis.

## Method reconciliation

| Pattern | Likely cause | Action |
|---|---|---|
| RAR allocation favours treatment, but time-adjusted effect shrinks | Secular trend was confounded with allocation | Make the time-adjusted pre-specified analysis primary and report drift sensitivity |
| BOIN and CRM choose different MTDs | Interval rule and model skeleton respond differently | Report both designs' OCs across calibrated dose-toxicity scenarios |
| A platform arm graduates but external Phase 3 fails | Selection and predictive-model optimism | Apply bias adjustment and evaluate calibration of the predictive threshold |
| CRP is claimed but simulated Type-I is high | Implementation does not preserve the formal conditional rejection probability | Re-derive the rule and validate it before use |

## Failure modes

### RAR drift bias

- Trigger: allocation changes while calendar effects, sites, or supportive care
  also change.
- Mechanism: time is confounded with treatment assignment.
- Symptom: the effect estimate is sensitive to calendar-time adjustment.
- Disposition: include pre-specified time covariates and proper analysis weights;
  report drift scenarios.

### RAR in a two-arm confirmatory trial

- Trigger: adaptive allocation is added for a presumed welfare advantage.
- Mechanism: extra bias and operational complexity provide little marginal
  benefit.
- Symptom: reviewers reject the justification or simulations reveal power loss.
- Disposition: use group-sequential stopping unless a setting-specific case is
  demonstrated.

### CRM skeleton miscalibration

- Trigger: a default skeleton is used without scenario calibration.
- Mechanism: the model's prior dose-toxicity spacing drives dose assignment.
- Symptom: systematic MTD-selection errors across plausible curves.
- Disposition: calibrate the skeleton or switch to a transparent interval design.

### Platform graduation bias

- Trigger: the selected arm's raw interim estimate is carried forward.
- Mechanism: selection favours positive noise.
- Symptom: subsequent studies show smaller effects.
- Disposition: report conditional/bias-adjusted estimates and predictive-
  probability calibration.
