# LCV (Latent Causal Variable)

## LCV (Latent Causal Variable)

LCV uses LDSC-merged genome-wide sumstats and reports gcp (genetic causality proportion) on [-1, 1]. It is a complement to, not a replacement for, MR; gcp ~ 0 with high LDSC rg implies pure genetic correlation without partial causation. LCV uses ALL genome-wide SNPs after LDSC-merging, not the MR instrument set.

```r
# RunLCV()'s own internal source("MomentFunctions.R") resolves against the R session's
# current working directory, not RunLCV.R's own location -- sourcing RunLCV.R without also
# cwd-ing into the clone's R/ folder crashes at call time with "cannot open file
# 'MomentFunctions.R'" (confirmed 2026-09-21, jean997/LCV HEAD). ldsc.intercept defaults
# to 1 (recommended for real data per LCV's own docs), which requires real n.1/n.2 -- pass
# your two GWAS sample sizes explicitly; the default n.1=n.2=1 silently mis-estimates.
setwd('LCV/R')  # assumes LCV was cloned as ./LCV, per this Skill's install step
source('RunLCV.R')
res_lcv <- RunLCV(ldscores$L2, x$Z, y$Z, n.1 = n_exposure, n.2 = n_outcome)
# res_lcv$gcp.pm (posterior mean gcp; there is no res_lcv$gcp field); res_lcv$pval.gcpzero.2tailed
```

**gcp interpretation:**

| gcp value | Interpretation |
|-----------|----------------|
| 0 | Pure genetic correlation; no partial causation |
| 0.5 | Partial causation; mixture |
| 0.6 | Partial causation; modestly causal direction |
| 1 | Fully causal in tested direction |
| Significant p_gcp != 0 | Directional evidence of (partial) causation |
