# Dispersion-mean trend does not fit

### Parametric dispersion-mean trend doesn't fit the cloud

**Trigger:** `plotDispEsts(dds)` shows the red parametric trend curve nowhere near the cloud of gene-wise (blue) and final (black) dispersion estimates.

**Mechanism:** Default `fitType='parametric'` assumes `dispersion ~ a/mean + b`. Fails when the experiment has very few samples per group with highly heterogeneous biology, many very-low-count genes pulling the trend, or a continuous covariate driving large variability.

**Symptom:** Curved trend that doesn't match the cloud; mismatch shows up most clearly in low-baseMean genes.

**Fix:** Refit with `DESeq(dds, fitType='local')` (local regression) or `fitType='mean'` (flat trend); compare `plotDispEsts` between fits and pick the one tracking the cloud. Falling back to `'mean'` is a sign the data is unusual; investigate before trusting results.
