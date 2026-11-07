# Background

Not needed to run an analysis. This explains choices the routes and scripts already make.

## Shrunken fold changes and p-values come from different models

`lfcShrink()` returns a fold change from a Bayesian posterior, but its `pvalue` is still the unshrunken Wald p-value from `results()`. The shrunken estimate is for ranking and plots; the p-value is for inference. Reporting both together is fine. Filtering on the shrunken fold change and then claiming false-discovery control for that filter is not: use `lfcThreshold=` in `results()` when a fold-change threshold must be part of the tested claim.

`results(dds)` with no `name=` or `contrast=` returns the last coefficient in `resultsNames(dds)`, which depends on factor level order and formula order. The scripts always pass an explicit contrast.

## Tests and estimators

| Test / estimator | What it tests | Limit |
|------------------|---------------|-------|
| Wald | One coefficient equals zero | Anti-conservative with many low-count outliers |
| LRT (`test='LRT'`, `reduced=`) | Joint effect of the dropped terms | Its fold-change column is the last coefficient only, not the joint effect |
| `lfcShrink(type='apeglm')` | Posterior fold change under a heavy-tailed prior | Needs `coef=`; cannot take `contrast=` |
| `lfcShrink(type='ashr')` | Posterior fold change under a unimodal prior | Accepts `contrast=`; slightly different inferential frame |
| `lfcShrink(type='normal')` | Posterior fold change under a zero-centred normal | Legacy; errors on designs with interactions |
| `lfcThreshold=` | Fold change exceeds a stated threshold | Conservative; use only for a pre-specified threshold |

## Small samples

At two or three samples per group every tool misses a substantial share of true changes. The result is still worth reporting, labelled exploratory.

## Filtering

Pre-filtering on total count is for speed and memory. Independent filtering happens later, inside `results()`, and sets low-count genes to `padj = NA` to keep power for the rest. Cook's distance flags single-sample outliers; it is not computed for continuous covariates, and at seven or more samples per group DESeq2 replaces outliers and refits.

## Pseudobulk

Summing counts per donor makes the donor the unit of replication. Testing cells as replicates treats thousands of correlated cells as independent and inflates false discoveries. Pseudobulk counts need none of the DESeq2 settings recommended for cells as columns (`test='LRT'`, `useT`, `minmu`, `minReplicatesForReplace=Inf`, glmGamPoi). DESeq2, edgeR quasi-likelihood and limma-voom agree closely on pseudobulk counts.
