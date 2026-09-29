# Regulatory Status, Review Questions, and References

Regulatory status changes. Verify current agency and ICH sources before using
these dated statements in a protocol, submission, or public claim.

## Dated regulatory landscape in the source skill

- FDA's *Adaptive Designs for Clinical Trials of Drugs and Biologics* became
  final on December 2, 2019. It discusses group-sequential designs, blinded and
  unblinded SSR, adaptive enrichment, and adaptive randomisation.
- FDA's oncology *Master Protocols* guidance became final in March 2022; 2018
  refers to the draft. Basket protocols study one drug across diseases,
  umbrella protocols study multiple drugs in one disease, and platform
  protocols allow arms to enter or exit over time.
- FDA separately published the broader *Master Protocols for Drug and
  Biological Product Development* revised draft in June 2026. Do not replace
  the still-final March 2022 oncology guidance with this broader draft; identify
  which document governs the statement being made.
- ICH endorsed E20 at Step 2 on June 25, 2025 and entered Step 3 public
  consultation. FDA issued the corresponding draft guidance in September 2025.
  It remained non-final when checked on September 28, 2026. Verify its present
  status before citing it.
- FDA's CDER Bayesian methodology guidance was a January 2026 draft in the
  source snapshot. It describes simulation-based Type-I calibration for
  Bayesian primary inference; verify current status before relying on it.
- Project Optimus (2021-2024) emphasises randomised dose comparison and the
  totality of dose/exposure evidence rather than an MTD-and-go strategy. FDA's
  dose-optimisation guidance became final in August 2024.
- BOIN received an FDA Fit-for-Purpose determination in December 2021. This is
  method qualification, not approval of a particular trial design.

## Quantities that must be contextualised

| Example quantity | Origin | How to use it |
|---|---|---|
| Promising-zone conditional power of about 30%-80% | Mehta-Pocock example | Treat as a pre-specified design choice and simulate alternatives |
| RAR most defensible with at least three arms | Robertson et al. consensus framing | Justify the setting; arm count alone is insufficient |
| OBF final nominal alpha near 0.024 of one-sided 0.025 | Typical four-look OBF-like design | Recalculate for the actual timing and spending function |
| I-SPY 2 graduation at predictive probability at least 85% for a 300-patient Phase 3 | I-SPY 2 example | Do not copy without calibrating the new platform |
| Pediatric power-prior discount around 0.3-0.6 | Example discussed in the source skill | Include in a sensitivity range, not as a default |
| Non-PH Schoenfeld underestimation of 20%-50% | Delayed-effect scenarios in the source skill | Simulate the trial-specific hazard trajectory |

## Anticipated reviewer questions

| Question | Evidence the response should contain |
|---|---|
| How is Type-I error controlled? | Exact closed-testing, spending, combination-test, or CRP construction plus null simulations and Monte Carlo error |
| Why these boundaries? | Clinical/operational rationale, information times, spending function, and final-analysis penalty |
| Was the SSR rule pre-specified? | Conditional-power zones, original weights, `n_max`, and the rule's implementation in the SAP |
| Is the IDMC firewall credible? | Access matrix, meeting/data-transfer SOP, and the exact message released to the sponsor |
| Why use RAR? | Multi-arm/rare-disease rationale, allocation floors, time-trend adjustment, analysis weights, and drift simulations |
| Why promising-zone rather than a CHW design? | Comparative OCs and a direct response to the Jennison-Turnbull critique |
| How is enrichment bias handled? | Selection rule, closed testing, raw and bias-adjusted estimators, and sensitivity analyses |
| Why BOIN rather than CRM? | Protocol usability plus comparative selection/overdose OCs under calibrated scenarios |

## Common citation and design errors

- Calling ICH E20 final without checking its current step.
- Citing the 2018 FDA Master Protocols draft as the final guidance.
- Describing Bauer-Kohne combination testing as obsolete despite its role in
  modern combination-test designs.
- Expecting an OBF design to stop easily at the first, low-information look.
- Claiming an unblinded sample-size increase preserves alpha without giving the
  weights or combination rule.
- Reporting selected-arm or enriched-population estimates without correction.
- Presenting simulated operating characteristics without random seeds,
  scenarios, replicate counts, or Monte Carlo uncertainty.

## Foundational references

- Babb J, Rogatko A, Zacks S. 1998. Cancer Phase I clinical trials: efficient
  dose escalation with overdose control. *Statistics in Medicine* 17:1103-1120.
- Bauer P, Kohne K. 1994. Evaluation of experiments with adaptive interim
  analyses. *Biometrics* 50:1029-1041.
- Berry DA. 2015. Commentary on Hey and Kimmelman. *Clinical Trials* 12:107-109.
- Cui L, Hung HMJ, Wang SJ. 1999. Modification of sample size in group
  sequential clinical trials. *Biometrics* 55:853-857.
- FDA. 2019. *Adaptive Designs for Clinical Trials of Drugs and Biologics*.
  Final Guidance.
- FDA. 2022. *Master Protocols: Efficient Clinical Trial Design Strategies to
  Expedite Development of Oncology Drugs and Biologics*. Final Guidance.
- FDA. 2026. *Master Protocols for Drug and Biological Product Development*.
  Revised Draft Guidance, June 2026.
- FDA. 2025. *E20 Adaptive Designs for Clinical Trials*. Draft Guidance,
  September 2025; based on ICH Step 2 endorsement in June 2025.
- FDA. 2026. *Use of Bayesian Methodology in Clinical Trials*. Draft Guidance
  in the source snapshot.
- Friede T, Kieser M. 2006. Sample size recalculation in internal pilot study
  designs. *Biometrical Journal* 48:537-555.
- Hey SP, Kimmelman J. 2015. Are outcome-adaptive allocation trials ethical?
  *Clinical Trials* 12:102-106.
- Jennison C, Turnbull BW. 2015. Adaptive sample size modification in clinical
  trials: start small then ask for more? *Statistics in Medicine*
  34(29):3793-3810.
- Lan KKG, DeMets DL. 1983. Discrete sequential boundaries for clinical trials.
  *Biometrika* 70:659-663.
- Liu S, Yuan Y. 2015. Bayesian optimal interval designs for Phase I clinical
  trials. *Journal of the Royal Statistical Society C* 64:507-523.
- Mehta CR, Pocock SJ. 2011. Adaptive increase in sample size when interim
  results are promising. *Statistics in Medicine* 30:3267-3284.
- Muller HH, Schafer H. 2001. Adaptive group sequential designs for clinical
  trials. *Biometrics* 57:886-891.
- O'Brien PC, Fleming TR. 1979. A multiple testing procedure for clinical
  trials. *Biometrics* 35:549-556.
- O'Quigley J, Pepe M, Fisher L. 1990. Continual reassessment method: a practical
  design for Phase I clinical trials in cancer. *Biometrics* 46:33-48.
- Pocock SJ. 1977. Group sequential methods in clinical trials. *Biometrika*
  64:191-199.
- Robertson DS, Lee KM, Lopez-Kolkovska BC, Villar SS. 2023.
  Response-adaptive randomization in clinical trials: from myths to practical
  considerations. *Statistical Science* 38:185-208.
- Wassmer G, Brannath W. 2016. *Group Sequential and Confirmatory Adaptive
  Designs in Clinical Trials*. Springer.
