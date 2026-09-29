---
name: bio-clinical-biostatistics-adaptive-designs
description: Designs adaptive clinical trials including group-sequential (O'Brien-Fleming, Pocock, Lan-DeMets spending), sample-size re-estimation (blinded Friede-Kieser, unblinded Cui-Hung-Wang, Mehta-Pocock promising zone), seamless Phase 2/3 with treatment-arm selection, population enrichment, and response-adaptive randomisation. Covers FDA 2019 Final Adaptive Designs Guidance, the March 2022 final oncology Master Protocols guidance, the June 2026 broader Master Protocols revised draft, and non-final ICH E20 guidance. Use when planning interim analyses, sample-size re-estimation, or master/platform-trial designs.
tool_type: r
primary_tool: rpact
goal_approach_exempt: true
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples were written for R with `rpact` 4.2+, `gsDesign` 3.6+,
`gsDesign2` 1.1+, `adaptr`, `simtrial`, `BOIN`, `dfcrm`, and `RBesT`.
Commercial alternatives include East/EastHorizon, ADDPLAN, and FACTS.

Before using a pattern, run `packageVersion("<pkg>")` and inspect the relevant
help page with `?<function>`. If an API differs, adapt the example to the
installed version rather than retrying unchanged. R is the regulatory de facto
standard; Python adaptive-design coverage is limited.

# Adaptive Clinical Trial Designs

Pre-specify every adaptation and preserve trial-wide Type-I error through a
valid group-sequential boundary, a combination test, closed testing, or the
Conditional Rejection Probability (CRP) principle. Do not treat an operational
rule as valid merely because it is called adaptive.

## Core Workflow

1. Define the estimand, endpoint, analysis population, sidedness, alpha, power,
   information scale, and non-proportional-hazards assumptions before choosing
   an adaptation.
2. Classify the proposed adaptation: early efficacy/futility stopping, blinded
   or unblinded sample-size re-estimation (SSR), treatment selection,
   enrichment, response-adaptive randomisation (RAR), or platform entry/exit.
3. Choose the Type-I control mechanism and pre-specify information times,
   boundaries or spending function, combination-test weights, decision zones,
   maximum sample size, multiplicity strategy, and any population-selection
   rule in the protocol and SAP.
4. Separate blinded and unblinded operations. For unblinded adaptations, define
   an independent data monitoring committee (IDMC) firewall and limit sponsor
   communication to the pre-specified operational decision.
5. Simulate operating characteristics under the null, design alternatives, and
   plausible misspecification. Report Type-I error, power, expected and maximum
   sample size, stopping probabilities, selection bias, and sensitivity to
   drift, delayed effects, nuisance parameters, and decision thresholds.
6. Produce a reproducible design object, executable simulation code, versioned
   software details, and SAP-ready language. Keep confirmatory claims tied to
   the simulated design actually implemented.

## Select the Design

| Scenario | Recommended approach | Why |
|---|---|---|
| Confirmatory early stopping | O'Brien-Fleming group-sequential design | Conservative early, small final-analysis penalty |
| Flexible actual look timing | Lan-DeMets alpha spending | Spending follows actual information fractions |
| Uncertain nuisance parameter | Blinded SSR (Friede-Kieser) | Preserves blinding and ordinarily avoids Type-I inflation |
| Increase sample size after an unblinded effect estimate | Inverse-normal/CHW combination design with a pre-specified promising zone | Uses original weights to preserve Type-I error |
| Seamless Phase 2/3 with arm selection | Bauer-Kohne combination test plus closed testing | Supports selection while preserving familywise error |
| Biomarker enrichment | Closed-test stage-wise enrichment | Controls testing across full and selected populations |
| Multi-arm rare-disease or biomarker platform | Bayesian platform/RAR with frequentist operating-characteristic calibration | Can support learning and allocation across at least three arms |
| Two-arm confirmatory trial considering RAR | Prefer group-sequential stopping | Avoids RAR drift and estimator bias with little welfare gain |
| Phase 1 dose finding | BOIN; compare against calibrated CRM/mTPI-2 OCs | BOIN has transparent pre-tabulated decisions |
| Phase 1b/2 dose optimisation | Randomised multi-dose design, such as BOIN-12 | Aligns with Project Optimus rather than MTD-and-go |
| Basket trial | EXNEX or robust MAP | Borrows while allowing a stratum to detach |
| Umbrella trial | Bayesian platform with shared control | Supports multiple therapies in one disease |

## Non-Negotiable Guardrails

- Do not naively increase sample size from an unblinded interim effect. Use
  pre-specified CHW/inverse-normal weights or another validated combination
  test and demonstrate Type-I control by simulation.
- Do not claim that a promising-zone rule has no cost. State its conditional
  power boundaries and maximum sample size and report simulated Type-I error;
  address the Jennison-Turnbull critique.
- Do not use RAR in a two-arm confirmatory setting by default. In multi-arm RAR,
  pre-specify calendar-time covariates, allocation floors, analysis weights,
  update cadence, and drift sensitivity analyses.
- Correct the selected-population and selected-arm estimates after enrichment
  or graduation; report raw and bias-adjusted estimates.
- Do not use a proportional-hazards Schoenfeld calculation for delayed-effect
  immunotherapy without a non-PH method such as Lakatos or simulation under the
  expected hazard-ratio trajectory.
- Never leak interim effect estimates through the IDMC firewall.
- Treat ICH E20 as non-final unless a current authoritative check establishes
  otherwise. Distinguish its June 2025 ICH Step 2 endorsement from FDA's
  September 2025 draft-publication date. Distinguish the March 2022 final
  oncology Master Protocols guidance from FDA's broader June 2026 revised
  draft and from the superseded 2018 oncology draft.

## Routed Guidance

- Read [group-sequential-ssr-enrichment.md](references/group-sequential-ssr-enrichment.md)
  for O'Brien-Fleming, Lan-DeMets, blinded/unblinded SSR, combination tests,
  CRP, and enrichment details.
- Read [platform-rar-dose-finding.md](references/platform-rar-dose-finding.md)
  for RAR, platform trials, BOIN/CRM, reconciliation, and method-specific
  failure modes.
- Read [regulatory-and-review.md](references/regulatory-and-review.md) before
  citing regulatory status, using quantitative thresholds, or preparing a
  reviewer response.
- Read [provenance.md](references/provenance.md) for source and license identity.
- Run or adapt [adaptive_designs.R](scripts/adaptive_designs.R) for the packaged
  executable group-sequential, blinded-SSR, BOIN, CRM, and MAP-prior examples.
  Promising-zone SSR, enrichment, platform graduation, and RAR are visibly
  documented-only scaffolds, not executable analyses. The MAP example is not
  an EXNEX implementation.
- Use [usage-guide.md](usage-guide.md) for prompt examples and the expected
  response shape.

## Required Deliverables

- Design rationale and explicit Type-I control mechanism.
- Protocol/SAP table of looks, information fractions, boundaries, spending,
  weights, adaptations, thresholds, sample-size caps, and analysis populations.
- Reproducible code plus package versions and random seeds.
- Simulation scenarios, Monte Carlo error, and operating-characteristic tables.
- IDMC access and communication plan for unblinded adaptations.
- Bias-adjusted estimates after selection or enrichment and a reconciliation
  when reasonable methods disagree.
- Clear separation of regulatory facts, design assumptions, and simulation
  results.

## Related Skills

- clinical-biostatistics/power-and-sample-size - sample size for adaptive designs
- clinical-biostatistics/multiplicity-graphical - closed testing
- clinical-biostatistics/bayesian-trials - Bayesian platforms and MAP priors
- clinical-biostatistics/trial-reporting - adaptive trial reporting
- clinical-biostatistics/survival-analysis - adaptive time-to-event designs
- experimental-design/sample-size - general sample-size methods
