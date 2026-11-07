# Choosing q


`q` is the latent-space (bottleneck) dimension of the fit; too low leaves confounders in the outlier signal, too high absorbs real signal, and it must stay far below the sample count. Estimate it per cohort instead of hard-coding it:

```r
# FRASER 2.6.1: estimateBestQ() (OHT); FRASER 2.2.0: optimHyperParams().
estimate_q <- if (exists('estimateBestQ', envir = asNamespace('FRASER'), inherits = FALSE)) FRASER::estimateBestQ else FRASER::optimHyperParams
fds <- estimate_q(fds, type = 'jaccard', plot = FALSE)
bestQ(fds, 'jaccard')
# Exhaustive injected-outlier grid is available with FRASER 2.6.1:
fds <- FRASER::estimateBestQ(fds, type = 'jaccard', useOHT = FALSE, q_param = c(2, 5, 10, 15), plot = FALSE)
plotEncDimSearch(fds, type = 'jaccard', plotType = 'auc')    # plotType = 'loss' also works
```

`plotEncDimSearch` without `plotType` shows the OHT singular-value plot and returns NULL (with a warning) after `useOHT = FALSE`. Measured on the 30-sample synthetic cohort (PCA): q=1-3 flagged 4/4 planted events, q=5-10 3/4 (the 60%-usage cryptic donor was missed), q=15 1/4; `estimateBestQ` chose q=1 (FRASER 2.6.1, OHT) and q=2 (2.2.0 grid), both 4/4. A fixed q=10 (the old default) is wrong for cohorts this small. OHT prints "Optimal latent space dimension is smaller than 2 ... set to 2" when the cohort is too small or degenerate (observed at n=4 in the synthetic check). For OUTRIDER, `estimateBestQ(ods)` (see its block) is the fast single-q estimate; `findEncodingDim` / `estimateBestQ(ods, useOHT = FALSE)` grids many q values and is slow.
