## fishpond/swish for Inferential-Uncertainty-Aware Testing

**Goal:** Test DTE while propagating quantification uncertainty from Salmon's Gibbs samples.

**Approach:** Run Salmon with `--numGibbsSamples 20`, import with tximeta, then use swish to average a non-parametric SAMseq-style test across inferential replicates.

```r
library(fishpond); library(tximeta)

# meta: data.frame(sample_id, condition); salmon_dir/<sample_id>/quant.sf from a run with --numGibbsSamples 20
coldata <- data.frame(names = meta$sample_id,
                      files = file.path(salmon_dir, meta$sample_id, 'quant.sf'),
                      condition = factor(meta$condition))   # swish needs a factor
se <- tximeta(coldata, skipMeta = TRUE)    # drop skipMeta once a linkedTxome matches your index
y <- scaleInfReps(se)
y <- labelKeep(y)
y <- y[mcols(y)$keep, ]

set.seed(1)
y <- swish(y, x = 'condition')
y <- computeInfRV(y)                        # adds mcols(y)$meanInfRV

dte_results <- as.data.frame(mcols(y))     # log2FC, pvalue, qvalue, meanInfRV
sig <- subset(dte_results, qvalue < 0.05)
```

Checked on real chrX Gibbs samples (fishpond 2.12.0, 2 v 2): 391 transcripts tested, 6 at q < 0.05 (top ENST00000380861, log2FC -1.45 with NumReads 75/63 vs 9/34); with a character `condition` swish stops with "is.factor(condition) is not TRUE". On a synthetic 6 v 6 set with planted DTE, swish called 33/35 planted up-transcripts (all log2FC > 0) with 29/848 false positives.

**`infRV`** (inferential relative variance) is a per-feature uncertainty diagnostic (`meanInfRV`); high-infRV transcripts are unreliable and can be filtered before testing (e.g. `y[mcols(y)$meanInfRV < 1, ]`, choose the cut from `summary(mcols(y)$meanInfRV)`; median 0.69 on the chrX set). Critical for genes with many similar isoforms (TTN, MAPT, NEFM) where Salmon's EM is uncertain.
