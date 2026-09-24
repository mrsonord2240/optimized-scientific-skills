# LFC shrinkage: apeglm, ashr, normal

## LFC Shrinkage -- Three Flavors, Three Failure Modes

**Goal:** Get a stable effect-size estimate appropriate for ranking, GSEA input, volcano-plot x-axis, and reporting.

**Approach:** Default to apeglm. Switch to ashr when arbitrary contrasts are needed. Never use normal for new analyses.

```r
res_apeglm <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')
res_ashr   <- lfcShrink(dds, contrast = c('condition','treated','control'), type = 'ashr')
```

Add `res = res` (from `results(..., alpha = ...)`) to either call to carry `padj` through unchanged; see the top of `SKILL.md`.

| Method | Prior | Accepts | Use when |
|--------|-------|---------|----------|
| apeglm | Cauchy (heavy-tailed) | `coef=` only | Default; preserves large effects, suppresses low-count noise |
| ashr | Unimodal scale-mixture | `coef=` or `contrast=` (incl. numeric) | Need contrast= for interaction sums or pairwise from `~ 0 + group` |
| normal | Zero-centered normal | `coef=` or `contrast=`; not interaction designs | Legacy only, to reproduce pre-1.16 shrunken LFCs (its p-values are still the unshrunken Wald ones) |

The apeglm-cannot-use-contrast footgun: if the question is "drug effect in KO" from `~ genotype * treatment`, apeglm cannot directly shrink that contrast. Workarounds: (a) rebuild as combined factor `~ 0 + group` and relevel so the desired comparison is a coefficient; (b) use ashr; (c) accept the unshrunken LFC for that one comparison.

`pvalue` does NOT change when shrinking (`lfcShrink()` preserves the Wald p-value from `results()`); `padj` does unless `res = res` is passed. See the Single Most Important Insight at the top.

### apeglm refuses arbitrary contrasts

**Trigger:** Question of the form "drug effect in KO genotype" from `~ genotype * treatment`; user tries `lfcShrink(dds, contrast=list(...), type='apeglm')`.

**Mechanism:** apeglm fits per-coefficient priors. A numeric or list contrast is a linear combination of coefficients, not a coefficient -- no prior to apply.

**Symptom:** Error: "type='apeglm' shrinkage only for use with 'coef'"

**Fix:** Rebuild design as `~ 0 + group` with `group = paste(genotype, treatment)`, relevel so the comparison is a coefficient, refit. Or use `type='ashr'` which accepts contrasts.
