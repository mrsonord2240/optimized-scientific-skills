# SIMEX correction for MR-Egger under NOME violation

Read when `I^2_GX < 0.9` and Egger must be reported. The runnable code is `scripts/simex_egger.R` (moved from SKILL.md "NOME violation invalidating Egger"); the trigger, mechanism and threshold stay there.

```bash
Rscript scripts/simex_egger.R --dat mr_out/harmonised.tsv --B 1000
```

`scripts/simex_egger.R` (simex) fits the weighted Egger `lm(beta.outcome ~ beta.exposure, weights = w, data = dat, x = TRUE, y = TRUE)` and runs `simex(SIMEXvariable = 'beta.exposure', measurement.error = dat$se.exposure, lambda = seq(0.5, 2, 0.5), fitting.method = 'quadratic', asymptotic = FALSE)`, then prints the SIMEX-corrected and the naive Egger slope. It precomputes the weights vector `w <- 1 / dat$se.outcome^2` rather than dividing a data-frame column in-formula: `simex()` refits the model on perturbed data and cannot re-evaluate `1 / se.outcome^2` against its own working frame, which has no `se.outcome` column, so the in-formula form crashes with "object 'se.outcome' not found" inside `simex()`'s refit (simex 1.8); passing a fully-qualified `data = dat` avoids it.
