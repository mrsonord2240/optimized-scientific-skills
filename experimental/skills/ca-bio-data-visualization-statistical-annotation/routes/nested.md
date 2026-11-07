# Nested observations

```bash
Rscript scripts/annotate_nested.R nested.csv nested.png holm
```

CSV requires nonmissing `group,subject_id,value`; each subject must belong to one group. The script fits a random-intercept model and annotates adjusted emmeans contrasts, saving the plot and a results CSV. Inspect model assumptions and convergence; it does not supply a full effect-size/uncertainty report. Repeated measures, factorial, survival and categorical outcomes need a design-specific model rather than this wrapper. Outputs may be overwritten.

Done when the model, biological replicate, comparison family, exact contrasts, effect magnitude and uncertainty have been checked.
