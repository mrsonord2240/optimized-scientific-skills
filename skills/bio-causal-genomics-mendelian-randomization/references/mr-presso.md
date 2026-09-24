# MR-PRESSO

Read when running the global / outlier / distortion tests or extracting outlier SNPs. The runnable code is `scripts/mr_presso_outliers.R` (moved from SKILL.md "MR-PRESSO Outlier Detection"); the prose below is verbatim.

## MR-PRESSO Outlier Detection

**Goal:** Detect horizontal-pleiotropy outliers, remove them, and test whether the corrected estimate differs from the uncorrected one (distortion test).

**Approach:** Three-step framework: global test (presence of pleiotropy), outlier test (per-SNP), distortion test (effect change after outlier removal).

```bash
Rscript scripts/mr_presso_outliers.R --dat mr_out/harmonised.tsv --nb 10000 --seed 42 --signif 0.05 --out mr_out/presso_outliers.txt
```

`scripts/mr_presso_outliers.R` fits `MRPRESSO::mr_presso()` (`OUTLIERtest`, `DISTORTIONtest`) on `dat[dat$mr_keep, ]` (`harmonise_data()` keeps dropped palindromes as `mr_keep = FALSE` rows; MR-PRESSO would fit them), seeds first (`set.seed(42)`: the global/outlier tests are Monte-Carlo, so the seed makes the p-value reproducible), prints the Global and Distortion tests, and writes the outlier SNPs to `--out`. `NbDistribution >= 10000` gives publication-grade p-value precision.

Outlier SNPs: `Outlier Test` is NULL unless the global test is significant. Its Pvalue is ALREADY Bonferroni-adjusted inside MRPRESSO (raw p x nrow(dat)) and may be a string such as "<3e-04", so do not divide the threshold by nrow(dat) again (that double correction flagged 0 outliers where MRPRESSO's own rule flagged 11, 9 of them planted; checked on MRPRESSO 1.0). MRPRESSO's own rule is adjusted P <= SignifThreshold, which the script applies (`p_adj <- as.numeric(sub('^<', '', ot$Pvalue))`, rows taken as `dat_p[rownames(ot)[which(p_adj <= 0.05)], 'SNP']`).

`NbDistribution` sets the cost and the floor: the outlier test needs `NbDistribution > nrow(dat) / SignifThreshold` (100 SNPs need > 2000), otherwise MRPRESSO warns "Outlier test unstable" and the outlier P is imprecise. 10000 draws on ~100 SNPs did not finish in 40 minutes on a shared CPU, so run it in the background, and use 3000-5000 while exploring.
