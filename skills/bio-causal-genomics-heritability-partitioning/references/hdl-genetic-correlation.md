# HDL genetic correlation (reference for heritability-partitioning SKILL.md)

## HDL Genetic Correlation (R)

**Goal:** Genetic correlation with ~60% lower variance than LDSC when samples are non-overlapping.

**Approach:** Reformat sumstats to HDL input -> download UKB reference panel eigen-decomposition -> run HDL.rg.

```r
library(HDL)

# Sumstats need columns: SNP, A1, A2, b (beta), se, N
gwas1 <- read.table('trait1.txt', header = TRUE, stringsAsFactors = FALSE)
gwas2 <- read.table('trait2.txt', header = TRUE, stringsAsFactors = FALSE)

LD.path <- 'UKB_array_SVD_eigen90_extraction'

rg_result <- HDL.rg(gwas1.df = gwas1, gwas2.df = gwas2, LD.path = LD.path,
                    Nref = 335265, output.file = 'hdl_rg.log',
                    eigen.cut = 'automatic')

# rg_result: rg, rg.se, p, h2_1, h2_2, gen.cov
```

HDL UKB reference (`UKB_array_SVD_eigen90_extraction`) requires non-overlapping samples; if either GWAS is from UKB, HDL is biased.

## Failure Mode: HDL bias with sample overlap

**Trigger:** Running HDL.rg on two GWAS that share > 5% of samples (e.g. two UKB traits, or UKB + FinnGen with overlapping recruitment).

**Mechanism:** HDL's eigen-decomposition likelihood treats the two traits as independent samples; sample-correlation in residuals biases the likelihood (typically inflates rg toward 1).

**Symptom:** HDL rg substantially different from cross-trait LDSC rg; HDL z-score very large compared to LDSC z; suspicious for high-correlation trait pairs.

**Fix:** Use cross-trait LDSC instead (intercept absorbs overlap). If both must be used, restrict HDL to unambiguously non-overlapping cohorts (e.g. UKB-only trait1 vs FinnGen-only trait2 with no shared individuals confirmed via individual ID exchange or IBD).

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| HDL rg = NA or numerical error | Sumstat format wrong; or eigen.cut too aggressive | Use eigen.cut='automatic'; check column names exactly |

## Install

```r
# HDL
remotes::install_github('zhenin/HDL/HDL')
# HDL UKB reference (downloads ~5 GB)
# https://github.com/zhenin/HDL/wiki/Reference-panels
```
