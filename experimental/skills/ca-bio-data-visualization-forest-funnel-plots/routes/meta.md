# Study-level forest and funnel

```bash
Rscript scripts/meta_plot.R studies.csv output-directory
```

CSV columns are `author`, `year`, `log_or`, and positive `log_or_se`; supply study-level log odds ratios and standard errors. The wrapper fits REML, uses HKSJ below five studies, reports I-squared only from five, and refuses pooling below three. It writes a ratio-scale forest with study and prediction intervals, a contour funnel, and a text summary; it runs Egger only at ten or more studies and trim-and-fill as sensitivity analysis. Existing named outputs are overwritten. Inspect the intervals and missing or invalid data before reporting. For k below three, plot individual estimates separately; that alternative is not implemented in this wrapper.

Done when effect scale, inputs, pooled interval, heterogeneity, prediction interval and limitations are checked; asymmetry alone does not prove publication bias.
