# Calibration demonstration

```bash
Rscript examples/targeted_quantification.R
```

This source example uses synthetic data, compares unweighted with 1/x² weighting, calculates back-calculated error and illustrates ion ratios; it writes to a temporary directory. It is not a user-data assay or a validated LLOQ calculation: choosing the lowest accuracy-passing point alone omits precision and full acceptance requirements. Adapt calibration and sample inputs, weighting, units, dilution, identity and validated range to the actual assay, with qualified review. R-squared alone does not establish low-end accuracy. The source's blank-area LOD formula mixes raw response with a ratio-calibration slope and must not be used without compatible units and a validated approach.

Done when the complete calibration, identity, accuracy, precision, range and reportability evidence have been reviewed for the intended use.
