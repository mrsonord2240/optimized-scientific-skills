# Adaptive Clinical Trial Designs - Usage Guide

## Overview

Design and stress-test pre-specified interim adaptations: group-sequential
stopping, blinded or unblinded sample-size re-estimation, treatment selection,
population enrichment, response-adaptive randomisation, platform entry/exit,
and Phase I dose finding.

## Prerequisites

R is the regulatory de facto standard. Install only the packages needed for the
chosen design, for example:

```r
install.packages(c(
    "rpact", "gsDesign", "gsDesign2", "adaptr", "simtrial",
    "BOIN", "dfcrm", "escalation", "trialr", "RBesT"
))
```

Commercial confirmatory-design tools include East/EastHorizon, ADDPLAN, and
FACTS. Verify package and guidance versions before use.

## Example prompts

### Group-sequential

> Design a one-sided 0.025 O'Brien-Fleming survival trial with looks at 33%,
> 67%, and 100% information, HR 0.65, 90% power, and explicit accrual,
> dropout, and minimum-follow-up assumptions. Return SAP-ready boundaries and
> simulation scenarios.

> Replace fixed look timing with a Lan-DeMets OBF-like spending function and
> show how actual information fractions alter cumulative alpha.

### Sample-size re-estimation

> Implement blinded variance SSR for a continuous endpoint at interim n=200,
> with an initial SD of 12 and n_max=300. State why the final test remains
> Type-I valid and simulate nuisance-parameter misspecification.

> Compare a Mehta-Pocock promising-zone design with an inverse-normal CHW
> design. Pre-specify conditional-power zones, original weights, n_max, and an
> IDMC firewall; report Type-I error and power with Monte Carlo uncertainty.

### Enrichment and platforms

> Build a closed-test enrichment design for full and biomarker-positive
> populations. Include the selection rule, multiplicity control, and a
> bias-adjusted selected-population estimate.

> Design a four-arm biomarker-stratified platform with allocation floors,
> response-adaptive updates, futility, graduation, calendar-time adjustment,
> and simulated secular drift.

### Phase I and dose optimisation

> Generate a BOIN escalation table for six doses, target DLT 30%, cohort size
> three, and maximum n=30. Compare MTD selection and overdose risk with a
> calibrated CRM over plausible dose-toxicity curves.

> Extend the selected Phase I doses into a Project Optimus-aligned randomised
> dose-comparison stage using safety, exposure, pharmacodynamics, and activity.

## Expected workflow

1. Define the estimand, endpoint, information scale, alpha, power, and plausible
   deviations from the working model.
2. Select a design and an explicit Type-I error-control mechanism.
3. Pre-specify every adaptation, threshold, weight, information time, and cap.
4. Define blinded/unblinded data access and the IDMC firewall.
5. Simulate null, alternative, and misspecified scenarios; include Monte Carlo
   error and selection/drift bias.
6. Deliver versioned executable code, operating-characteristic tables, and
   protocol/SAP-ready language.

## Tips

- Group-sequential stopping is generally easier to defend than RAR in a two-arm
  confirmatory trial.
- Blinded and unblinded SSR are different designs; never swap them after seeing
  interim results.
- Lan-DeMets spending gives timing flexibility, but the actual information
  fractions still drive the boundaries.
- Treat promising-zone thresholds and platform graduation probabilities as
  calibrated design choices, not reusable constants.
- Correct selected-arm and selected-population estimates for selection bias.
- Verify current ICH E20 and FDA guidance status before making regulatory claims.

## Routed resources

- [SKILL.md](SKILL.md) - mandatory workflow, design choice, and guardrails
- [Group-sequential, SSR, and enrichment](references/group-sequential-ssr-enrichment.md)
- [RAR, platforms, and dose finding](references/platform-rar-dose-finding.md)
- [Regulatory status and reviewer preparation](references/regulatory-and-review.md)
- [Runnable R examples](scripts/adaptive_designs.R) - executable
  group-sequential, blinded-SSR, BOIN, CRM, and MAP demonstrations plus clearly
  labelled documented-only scaffolds
- [Focused regression harness](scripts/test_adaptive_designs.R) - sources the
  combined runner and checks the repaired SSR, BOIN, CRM, MAP, and
  documented-only contracts
