# Group-Sequential, Sample-Size Re-estimation, and Enrichment Designs

Use this reference after selecting one of these adaptation families. The values
below are examples, not universal defaults; pre-specify and calibrate the actual
trial design.

## Group-sequential designs

O'Brien-Fleming (OBF) boundaries are very conservative at early looks and near
nominal at the final analysis. Pocock boundaries use a nearly constant nominal
alpha and therefore impose a larger final-analysis penalty. A Lan-DeMets
spending function preserves the intended spending shape when analyses occur at
different information fractions than planned.

```r
library(gsDesign)

design_obf <- gsDesign(
    k = 4,
    test.type = 1,
    alpha = 0.025,
    beta = 0.10,
    sfu = sfLDOF,
    timing = c(0.25, 0.50, 0.75, 1.0)
)
print(design_obf)
plot(design_obf)
```

For a time-to-event design, make the accrual, dropout, control hazard, hazard
ratio, total duration, and minimum follow-up assumptions explicit:

```r
n_gs <- gsSurv(
    k = 3,
    test.type = 2,
    alpha = 0.025,
    beta = 0.10,
    sfu = sfLDOF,
    lambdaC = 0.04,
    hr = 0.70,
    eta = 0.005,
    T = 24,
    minfup = 12
)
```

Do not use this proportional-hazards calculation unchanged for a delayed
immunotherapy effect. Use Lakatos or trial simulation under the expected
time-varying hazard ratio.

## Blinded SSR

Friede-Kieser blinded SSR re-estimates a nuisance parameter such as variance,
control event rate, or overall event rate without using the treatment-effect
estimate. A pre-specified nuisance-only rule can preserve Type-I error when the
recalculation is genuinely blinded, no treatment-effect information leaks into
the decision, and the valid final test is unchanged. This is conditional, not
an unconditional guarantee. Define the recalculation time, estimator, minimum
and maximum sample sizes, operational blinding, and final test.

```r
library(rpact)

design_blinded_ssr <- getDesignGroupSequential(
    kMax = 2,
    alpha = 0.025,
    beta = 0.20,
    sided = 1,
    informationRates = c(0.5, 1)
)
```

## Unblinded SSR and promising zones

Naively increasing sample size after seeing an interim treatment effect inflates
Type-I error. Cui-Hung-Wang preserves alpha with pre-specified weights based on
the original design:

```text
Z_weighted = w_1 * Z_1 + w_2 * Z_2_residual
```

Use an inverse-normal combination design or another validated implementation:

```r
design_unblinded_ssr <- getDesignInverseNormal(
    kMax = 2,
    alpha = 0.025,
    beta = 0.20,
    sided = 1,
    informationRates = c(0.5, 1),
    typeOfDesign = "WT",
    deltaWT = 0.25
)
```

A Mehta-Pocock rule often divides conditional power into:

- unfavourable, below about 30%: stop or continue without increasing;
- promising, about 30% to 80%: increase up to a pre-specified maximum;
- favourable, above about 80%: continue without increasing.

These thresholds are design choices, not defaults. The promising-zone
construction has been described as producing only small unconditional Type-I
inflation, but Jennison and Turnbull criticised this as hidden or "stealth"
alpha. Report the actual simulated error rate and compare the design with a
CHW-weighted group-sequential alternative.

The IDMC receives the unblinded effect estimate. The sponsor receives only the
pre-specified decision, such as "increase" or "no increase". Document the
firewall and any information that can be inferred from that decision.

The packaged `scripts/adaptive_designs.R` does not implement this full rule and
labels its promising-zone section `DOCUMENTED_ONLY`. Before making it
executable, supply conditional power, the adaptation function, `n_max`, the
original combination-test weights, IDMC output, and calibrated null and
alternative simulations with Monte Carlo uncertainty.

## Combination tests and CRP

Bauer-Kohne combines independent stagewise p-values, commonly through Fisher's
product. Inverse-normal combination tests use pre-specified stage weights.
Muller-Schafer's Conditional Rejection Probability principle permits a design
change when the null conditional rejection probability is preserved.

```r
design_fisher <- getDesignFisher(
    kMax = 3,
    alpha = 0.025,
    sided = 1
)

design_inverse_normal <- getDesignInverseNormal(
    kMax = 3,
    alpha = 0.025,
    informationRates = c(0.33, 0.67, 1.0)
)
```

Do not invoke the CRP principle by name without verifying the implemented
conditional rejection function. An ad hoc rule can still fail simulation.

## Adaptive enrichment

Define the full population and candidate subgroup before the trial. Use closed,
stagewise testing to control familywise error when a negative subgroup may be
dropped and the trial re-powered in a biomarker-positive subgroup. Simulate the
selection rule under null and heterogeneous alternatives.

Selection inflates the observed effect in the chosen population. Pre-specify a
bias-adjusted estimator, such as conditional maximum likelihood or a calibrated
hierarchical Bayesian approach, and report both raw and adjusted estimates.

## Reconciliation checks

| Disagreement | Likely cause | Required response |
|---|---|---|
| Blinded and unblinded SSR propose different sample sizes | One uses only a nuisance parameter; the other uses the treatment effect | Keep the pre-specified method; do not switch after inspection |
| A naive group-sequential result and CHW final test disagree | Different information weights or test definitions | Apply the pre-specified boundary and weights |
| Promising-zone and CHW designs propose different increases | Different calibration and efficiency assumptions | Compare operating characteristics under the same scenarios |
| Enriched-population replication is smaller | Interim selection bias | Report bias-adjusted estimates and replication sensitivity |

## Failure signatures

- Type-I error above target after an unblinded increase: weights or combination
  logic were absent or implemented incorrectly.
- Conditional power boundaries or `n_max` absent from the SAP: the SSR rule was
  not fully pre-specified.
- Treatment effect leaked beyond the IDMC: the unblinded adaptation is
  operationally compromised.
- Selected-population estimate reported without adjustment: the enrichment
  result is vulnerable to winner's-curse bias.
