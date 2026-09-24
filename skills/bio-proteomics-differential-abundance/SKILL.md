---
name: bio-proteomics-differential-abundance
category: Data Analysis
description: Tests for differentially abundant proteins between conditions with limma/DEqMS empirical-Bayes moderation, proDA/msqrob2/MSstats missingness modeling, and Python Welch+BH alternatives. Frames missing values as left-censored MNAR (model, do not impute), makes variance moderation the load-bearing step at n=3-5, and prefers feature/peptide-level testing. Use when identifying proteins with significant abundance changes between experimental groups. Summarization and normalization mechanics are proteomics/quantification; volcano and MA plots are data-visualization/volcano-and-ma-plots; pathway enrichment of the hit list is pathway-analysis/go-enrichment.
tool_type: mixed
primary_tool: limma
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: limma 3.62.2, DEqMS 1.24.0, proDA 1.20.0, msqrob2 1.14.1, QFeatures 1.16.0, MsCoreUtils 1.18.0, MSstats 4.14.2, ashr 2.2.63 (R 4.4.3 / Bioconductor 3.20), pandas 3.0.5, scipy 1.18.1, statsmodels 0.15.0 (checked 2026-09-15; feature-level material in `references/feature_level.md` rechecked 2026-09-21)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Differential Protein Abundance -- Moderated Testing on a Log-Intensity Matrix with Honest Missingness

**"Find differentially abundant proteins between my conditions"** -> Moderated statistical testing on a normalized log-intensity matrix, carrying missingness in the likelihood instead of filling it in -- because the missing values are low BECAUSE the protein is low, and imputing them replaces a detection limit with numbers whose fold change is set by the imputation constant and whose p-value carries no FDR guarantee.
- R: `limma::eBayes(fit, trend=TRUE, robust=TRUE)` for empirical-Bayes moderated t-tests (the protein-level workhorse)
- R: `DEqMS::spectraCounteBayes()` when PSM/peptide counts are available (preferred over limma-trend when quant depth varies)
- R: `proDA::test_diff()` when missing values are extensive (model the dropout, no imputation)
- R: `msqrob2::msqrob()` on a `QFeatures` object, or `MSstats::groupComparison()`, to test from the peptide/precursor table instead of a protein matrix
- Python: `scipy.stats.ttest_ind(equal_var=False)` + `statsmodels` BH (large n only; no moderation)

Scope: this skill owns the statistical TEST on a log2-transformed, normalized matrix (or a peptide table) -- design/contrast construction, variance moderation, missingness handling, multiple-testing correction, minimum-fold-change testing, and fold-change shrinkage. Peptide-to-protein summarization and normalization mechanics -> proteomics/quantification. Volcano/MA plots -> data-visualization/volcano-and-ma-plots. Enrichment of the hit list -> pathway-analysis/go-enrichment. OUT OF SCOPE: how MaxLFQ/TMP/IRS produce the matrix (quantification); how to draw a volcano (data-visualization); single-sample (n=1) comparisons, and diagnosis or treatment decisions for an individual patient.

## The Single Most Important Modern Insight -- Model the Missingness, Moderate the Variance, Test at the Feature Level

1. **Missing values in label-free MS are left-censored MNAR -- missing BECAUSE the intensity is low -- and imputing them (especially Perseus/MaxQuant downshift) turns a detection limit into numbers the test treats as measurements. Downshift draws each missing value from a Gaussian with mean = mu - 1.8*sigma and SD = 0.3*sigma, where mu and sigma describe ALL proteins in that run (sigma is not the replicate SD). For an on/off protein (seen in all of group A, missing in all of B) the fold change is therefore set by the imputation constant and the protein's own intensity, not by biology -- the volcano-plot "anchor/wing" artifact (on/off proteins on a streak whose logFC tracks average intensity). The p-value carries no FDR guarantee: whether 0.3*sigma is narrower or wider than the real replicate SD decides whether imputation inflates false positives or only costs power. The honest statement is "undetected in group B", not "20x lower, p=1e-6". The correct approach is to MODEL the dropout in the likelihood -- proDA (probabilistic dropout), msqrob2, MSstats-AFT -- NOT fill it (Lazar 2016; Ahlmann-Eltze & Anders 2019).
2. **At the n=3-5 replicates proteomics actually uses, per-protein variance has only 2-4 residual df and is unusable raw -- variance moderation is the load-bearing element, not optional.** limma borrows a prior d0 across all proteins so a 4-replicate design tests on ~10 df instead of 6; `trend=TRUE` makes the prior a function of mean intensity (effectively mandatory for label-free, where a single global prior mis-calibrates FDR across the abundance range); `robust=TRUE` Winsorizes outlier variances (Phipson 2016). DEqMS makes the prior a function of PSM/peptide count and generally outperforms limma-trend when quantification depth varies across proteins (Zhu 2020).
3. **Feature/peptide-level modeling keeps information that summarize-then-test throws away -- but it is not free.** Summarizing first (one number per protein per run) discards the within-protein between-peptide variance and the correct degrees of freedom: 12 consistent peptides deserve a smaller SE than 12 disagreeing ones, but after summarization both look equally certain, and a protein with 30 observations looks as informative as one with 3. msqrob2/MSstats keep every peptide as a degree of freedom (Goeminne 2016; Sticker 2020; Choi 2014). Two calibrations from running all three on the same 4 v 4 peptide table: on clean, balanced data the feature-level fit and limma on the robustly summarized matrix returned the SAME 80 calls -- the gain is real when peptides within a protein disagree or coverage is unbalanced, not automatic. And the smaller SEs cut both ways: per-run median normalization of the peptide table, which a protein-level test mostly shrugs off, produced 31 false positives in 111 calls at the feature level (it adds a between-condition offset AND shrinks the residual), so a feature-level result is only readable after the centring checks in `references/feature_level.md`.

## Tool Taxonomy

| Tool / method | Citation | Mechanism / role | When |
|---------------|----------|------------------|------|
| limma | Ritchie 2015; Phipson 2016 | EB moderated t; posterior variance blends a prior d0 with the per-protein estimate; `trend` ties the prior to mean intensity, `robust` Winsorizes outliers | protein-level summaries, small n, the default workhorse |
| DEqMS | Zhu 2020 | prior variance = loess of log-variance vs log2(count); precision follows quantification DEPTH not just intensity | TMT (count=PSM) and label-free DDA (count=peptide); quant depth varies; preferred over limma-trend |
| proDA | Ahlmann-Eltze & Anders 2019 (preprint) | probabilistic dropout: missing = left-censored, integrated under a per-sample sigmoid dropout curve; EB on location and variance; no imputation | label-free DDA with many MNAR missing values, small n, proteins absent in one group |
| msqrob2 | Sticker 2020; Goeminne 2016 | peptide-level robust ridge on a `QFeatures` object: Huber M-estimation downweights outlier peptides, ridge shrinks effects from few observations (>2 mean-model parameters only), EB variance moderation | label-free DDA, outlier-peptide / unbalanced-coverage risk; best FDR in hard spike-in regimes |
| MSstats | Choi 2014 | feature-level linear mixed model (group fixed + feature + run/subject random); AFT censored imputation only inside `dataProcess(MBimpute = TRUE)` (`groupComparison` has no censoring argument) | SRM/PRM/DIA, technical replicates, nested/repeated-measures, labeled designs |
| Welch t-test + BH | -- | per-protein two-sample t with `equal_var=False` + Benjamini-Hochberg | large n (>10/group), Python-only; no moderation, unusable at n=3-5 |
| ashr | Stephens 2017 | mixture prior with a point mass at zero; posterior means shrink uncertain effects toward zero | recovering "which proteins truly changed and by how much" (not for GSEA ranking) |
| volcano / MA plot | -- | (route OUT) | visualization -> data-visualization/volcano-and-ma-plots |
| enrichment of hits | -- | (route OUT) | functional interpretation -> pathway-analysis/go-enrichment |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Small n (3-5/group), protein-level summary matrix | limma `eBayes(trend=TRUE, robust=TRUE)` | EB borrows variance across proteins; the trend calibrates FDR across abundance |
| PSM/peptide counts available (TMT or label-free DDA) | DEqMS `spectraCounteBayes` | prior keyed on quant depth removes single-PSM false positives limma admits |
| Label-free with many MNAR missing values, on/off proteins | proDA `test_diff` | models the censored dropout; never imputes. It is honest but underpowered for on/off proteins at n=3-5 (0 of 11 called at n=4, best adj_pval 0.17): report those as a separate undetected list, not as non-significant |
| A peptide/precursor table is available, not just a protein matrix | msqrob2 (`QFeatures` + `msqrob`) or MSstats; workflows in `references/feature_level.md` | the feature table is the only place the between-peptide spread still exists; both keep it |
| Outlier-peptide risk, unbalanced peptide coverage | msqrob2 `robustSummary` + `msqrob(robust = TRUE)` | Huber M-estimation downweights the bad peptide instead of letting it move the protein mean. On a clean, balanced 4 v 4 it buys nothing: msqrob2 and limma on the same summarized matrix returned the identical 80 calls. Escalate for the disagreement, not by default |
| Three or more groups, or group plus covariates, at the peptide level | msqrob2 `msqrobAggregate(..., ridge = TRUE)` | ridge shrinks effects estimated from few observations -- but it is REFUSED on a two-group design (needs >2 mean-model parameters); use `ridge = FALSE` there |
| Technical replicates, nested/repeated-measures, labeled (SRM/PRM/DIA) | MSstats (feature-level mixed model) | random effects capture run/subject structure summarize-then-test discards |
| Any feature-level route | run the two centring checks in `references/feature_level.md` before reading the table | per-run median normalization of the peptide table (offset -0.18 plus a halved SE) gave 20.8-27.9% realized FDR at nominal 5% |
| Batch present | batch as a covariate in the design (`~ batch + condition`) | `removeBatchEffect` is visualization-only; never feed its output to `lmFit` |
| Minimum biologically meaningful fold change | `treat()` + `topTreat()` (or SAM s0) | tests |log2FC|>c against the moderated null; a post-hoc FC+significance double filter inflates FDR |
| Large n (>10/group), Python-only | Welch t-test + BH | variance estimates reliable; no moderation needed at large n |

Default when uncertain: protein-level summary matrix at n=3-5 -> limma `eBayes(trend=TRUE, robust=TRUE)`; if PSM/peptide counts exist, escalate to DEqMS; if missingness is extensive and intensity-dependent, escalate to proDA.

## limma Workflow (R)

**Goal:** Identify differentially abundant proteins using moderated statistics that borrow information across all proteins.

**Approach:** Filter to proteins with enough valid values per group (rows with no or too few values give `NA` averages that stop `eBayes(trend = TRUE)`), build the design (batch as a covariate when present), fit the linear model, drop rows the design cannot estimate, then contrast, apply EB moderation with the intensity trend and robust fitting, and extract BH-corrected results. The per-condition count is NOT sufficient in a paired or blocked design (donor, subject or batch in the model): rows whose observed values fall in different blocks in the two conditions leave zero residual df or a partly `NA` coefficient vector, and limma tests them anyway on a contrast that is partly a block difference. Report proteins removed by either filter (e.g. undetected in one group) separately. Never feed `removeBatchEffect` output to `lmFit`. Runnable end to end on simulated data: `examples/limma_analysis.R`.

Runnable script: `scripts/limma_de.R` holds the valid-value filter, the design with batch as a covariate, the estimability filter, `eBayes(trend = TRUE, robust = TRUE)` and the BH table. `matrix.csv` is log2, proteins by samples; `samples.csv` has `sample`, `condition` and optionally `batch`; the contrast names are levels of `condition`. Its output lists the rows dropped by each filter. Report them.

```bash
Rscript scripts/limma_de.R matrix.csv samples.csv Treatment-Control results.csv [min_valid=2]
```

To keep the fit for `treat()`, DEqMS and ashr below, `source('scripts/limma_de.R')` and call `out <- run_limma_de(protein_matrix, sample_info, 'Treatment-Control')`, then `fit2 <- out$fit2` and `results <- out$results`.

### Minimum-fold-change testing

**Goal:** Call proteins whose effect exceeds a biologically meaningful threshold, not merely differ from zero.

**Approach:** Use `treat()` against the moderated null and read `topTreat()`. `treat()` re-estimates the prior and its `trend`/`robust` default to FALSE, so pass both again. Keep the result in its own object so `fit2` stays the eBayes fit for DEqMS and ashr. NEVER `topTable(lfc=...)` nor a post-hoc volcano double filter (`abs(logFC) > 1 & adj.P.Val < 0.05`): the BH guarantee then refers to FC = 0, not to the FC threshold, and conditioning on both FC and p can select high-variance nulls. How much FDR inflates depends on the regime -- above 50% for top-ranked lists with many small effects (Ebrahimpoor & Goeman 2021), far less when true effects are large.

```r
LFC_THRESHOLD <- log2(1.2)  # 1.2-fold floor; treat tests against this null, no double-filter FDR inflation
fit_treat <- treat(fit2, lfc = LFC_THRESHOLD, trend = TRUE, robust = TRUE)  # treat's trend/robust default to FALSE
results <- topTreat(fit_treat, coef = 1, number = Inf)  # topTreat omits the B column
```

## DEqMS Workflow (R)

**Goal:** Improve on limma by tying each protein's prior variance to its quantification depth -- proteins measured by more PSMs/peptides are more precise.

**Approach:** Run limma through `eBayes`, attach the count vector, then apply DEqMS's count-aware EB. Use PSM count for TMT (quant at MS2) and peptide count for label-free DDA; for multi-batch TMT use the MINIMUM count across batches (the bottleneck batch sets precision).

```r
library(DEqMS)

# fit2 is the limma fit through eBayes (above), after the valid-value filter
stopifnot(all(fit2$df.residual > 0))  # NA-sigma rows make spectraCounteBayes recycle loess predictions onto the wrong proteins
fit2$count <- psm_count_per_protein[rownames(fit2$coefficients)]  # PSM for TMT, peptide for LFQ; min across batches
fit3 <- spectraCounteBayes(fit2)

results <- outputResult(fit3, coef_col = 1)
# adds sca.t, sca.P.Value, sca.adj.pval (the count-adjusted statistics; use these, not the limma columns)
```

## proDA Workflow (R)

**Goal:** Test proteins with extensive MNAR missingness, including on/off proteins, without imputing a single value.

**Approach:** Fit the probabilistic-dropout model directly on the log-intensity matrix; missing values contribute as left-censored observations under a per-sample dropout curve. Then test the contrast against zero. proDA keeps on/off proteins in the model honestly but rarely reaches significance for them at n=3-5 -- a censored observation carries less information than a measured one -- so list them as "undetected in group X" rather than reading their non-significance as evidence of no change.

```r
library(proDA)

fit <- proDA(protein_matrix, design = ~condition + batch, col_data = sample_info,
             reference_level = 'Control')
result_names(fit)  # Intercept, conditionTreatment, batch...: test a coefficient name
results <- test_diff(fit, 'conditionTreatment')
# columns: name, pval, adj_pval, diff (log2FC), t_statistic, se
# do not report diff for proteins with no observed value in a group: the location prior sets it (sign can be wrong)
```

## Feature-Level Workflows (msqrob2, MSstats)

Read `references/feature_level.md` before testing from a peptide/precursor table. It holds the runnable msqrob2 blocks (protein-summary and peptide-level `msqrobAggregate`), the MSstats block, the two centring checks, and their Common Errors. Rules that hold whichever route you take:
- Report proteins with no observation in a condition as undetected: msqrob2 returns `adjPval = NA` for them, MSstats an infinite `log2FC` with `issue == 'oneConditionMissing'`. Neither is a fold change or "not significant".
- Do not read a feature-level table until `median(log2FC)` and the residual-SD check have been run; per-run median normalization on the peptide table is the failure they catch. End-to-end example: `examples/msqrob2_peptide_level.R`.

## Python Workflow

**Goal:** Run the full pipeline in Python when no R is available and n is large enough that moderation is unnecessary.

**Approach:** Log2-transform, median-normalize, run per-protein Welch t-tests, apply Benjamini-Hochberg. Return the untestable proteins (fewer than 2 values in a group) alongside the results, as the limma section requires -- they are not "not significant". This has NO variance moderation and should not be used at n=3-5 -- escalate to limma/DEqMS for small n. Runnable script with the `preprocess()` and `differential_abundance()` functions: `examples/differential_abundance.py` (run it as `python examples/differential_abundance.py`; import the two functions for your own matrix).

## Fold-Change Reporting

**Goal:** Hand the right effect estimate to the right consumer.

**Approach:** Report the RAW fold change (the best unbiased point estimate) for GSEA/pathway ranking and meta-analysis -- those need the full continuous distribution or FC+SE pairs. Apply shrinkage (ashr) only when recovering "which proteins truly changed and by how much"; it fits a mixture prior with a point mass at zero and shrinks uncertain effects smoothly toward zero. This is preferred over hard-thresholding (zeroing FCs at padj 0.05), which creates an arbitrary step function. For a reporting table give the raw FC with its adjusted p and confidence interval; for meta-analysis pool raw FCs with their standard errors and shrink afterwards, not before. No mature Python ashr equivalent exists.

```r
library(ashr)

ok <- !is.na(fit2$coefficients[, 1]) & !is.na(fit2$s2.post)  # ash() returns PosteriorMean 0 / prior lfsr for NA rows
se <- sqrt(fit2$s2.post[ok]) * fit2$stdev.unscaled[ok, 1]
shrunk <- ash(fit2$coefficients[ok, 1], se, mixcompdist = 'normal')
shrunken_fc <- shrunk$result$PosteriorMean  # report alongside raw logFC, not as a replacement for GSEA
lfsr <- shrunk$result$lfsr
```

## Per-Method Failure Modes

### Downshift / any imputation feeding a variance-based test
**Trigger:** Perseus/MaxQuant downshift (or MinDet/MinProb/QRILC) fills NAs, then limma/t-test runs on the filled matrix.
**Mechanism:** Insight 1 -- the on/off fold change is fixed by the downshift constant and the imputed SD is unrelated to the real replicate SD, so the test is anticonservative or merely underpowered depending on the data.
**Symptom:** Volcano "anchor/wing" -- on/off proteins on a streak whose logFC tracks average intensity (rigid near-vertical streaks of pinned points are imputation artifacts, not biology); p-values and FDR without a guarantee.
**Fix:** Model the missingness instead (proDA / msqrob2 / MSstats-AFT); report on/off proteins as "undetected in group X".

### kNN imputation on left-censored data
**Trigger:** kNN/mean imputation applied to label-free data with MNAR dropout.
**Mechanism:** Mean-reverting -- pulls a truly-low (missing because low) value UP toward the mean.
**Symptom:** Real down-regulation is compressed; down hits weakened or lost.
**Fix:** Only valid under MCAR/MAR; for MNAR model the dropout. Under uncertainty Lazar 2016 shows the milder MCAR error beats MNAR-imputers slamming random highs to the floor.

### removeBatchEffect before testing
**Trigger:** `removeBatchEffect()` output fed to `lmFit`.
**Mechanism:** Subtracts the fitted batch component with no uncertainty propagation -> understated residual variance, inflated EB df; if batch is confounded with biology it deletes real signal.
**Symptom:** Anticonservative p-values; lost true effects when cases/controls split by batch.
**Fix:** Include batch as a covariate in the SAME model (`~ batch + condition`); use `removeBatchEffect` only for PCA/visualization.

### eBayes(trend=FALSE) on intensity data
**Trigger:** Plain `eBayes` (trend off) on a log-intensity matrix.
**Mechanism:** A single global prior over-shrinks high-abundance and under-shrinks low-abundance proteins.
**Symptom:** Mis-calibrated FDR across the abundance range.
**Fix:** `eBayes(trend = TRUE, robust = TRUE)`; escalate to DEqMS when quant depth varies.

### Wrong DEqMS count column
**Trigger:** Razor+unique counts vs MS2-level PSMs, or total-across-batches vs minimum-across-batches.
**Mechanism:** The variance-vs-count prior is fit on the wrong precision proxy.
**Symptom:** Mis-ranked proteins; the count moderation helps the wrong ones.
**Fix:** PSM count for TMT, peptide count for label-free; minimum count across batches for multi-batch TMT.

### proDA on MCAR missingness
**Trigger:** proDA applied where dropout is random (e.g. a TMT channel lost at random), not detection-limited.
**Mechanism:** The left-censored dropout model is mis-specified.
**Symptom:** Biased estimates; the model fits a dropout curve that does not exist.
**Fix:** proDA needs intensity-dependent missingness; for MCAR use limma/DEqMS on the observed values.

Feature-level failure modes (per-run median normalization of a peptide table) are in `references/feature_level.md`. The FC + significance double filter is covered under Minimum-fold-change testing.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| min-FC floor log2(1.2) (1.2-fold) via treat() | -- | example floor; common alternatives 1.5-fold (~0.58) or 2-fold (1.0); set by biology, tested against the moderated null |
| BH adjusted p < 0.05 | Benjamini-Hochberg | controls FDR over the WHOLE rejection set, not subsets carved out afterward |
| feature-level centring trip-wires (`OFFSET_MAX` 0.05, `SD_RATIO_MIN` 0.85) | measured, one 4 v 4 set | empirical, not distributional bounds; in `references/feature_level.md` |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `eBayes`: `prior.weights contain NA values` | rows with no valid value (MaxQuant all-zero LFQ rows) or too few per group | valid-value filter before `lmFit` |
| `lmFit`: `Partial NA coefficients for N probe(s)`; `stopifnot(all(fit2$df.residual > 0))` fires | a paired/blocked design (donor, subject, batch) where a row's values sit in different blocks per condition; the per-condition valid-value filter does not catch this | `fit <- fit[fit$df.residual > 0 & rowSums(is.na(fit$coefficients)) == 0, ]` after `lmFit` |
| DEqMS warning `longer object length is not a multiple of shorter object length` | rows with NA sigma / zero residual df reach `spectraCounteBayes` | filter before `lmFit`; require `fit2$df.residual > 0` |
| proDA: `object 'conditionControl' not found` | intercept design with `reference_level`: coefficients are `Intercept`, `conditionTreatment` | `test_diff(fit, 'conditionTreatment')` |
| `makeContrasts`: `object 'DrugB' not found` | design rename hard-coded to two groups | `colnames(design)[seq_len(nlevels(cond))] <- levels(cond)` |
| Python `ValueError: No protein has >= 2 non-missing values` | a group has one sample, or no protein is observed twice per group | not a two-sample test; do not run DA on n=1 |
| DEqMS columns missing | forgot `fit$count` or read limma columns | set `fit$count`, run `spectraCounteBayes`, read `sca.adj.pval` from `outputResult` |

Feature-level errors (QFeatures, msqrob2, MSstats, MaxQuant flag columns) are in `references/feature_level.md`.

## References

- Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. 2015. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res* 43(7):e47.
- Phipson B, Lee S, Majewski IJ, Alexander WS, Smyth GK. 2016. Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression. *Ann Appl Stat* 10(2):946-963.
- Zhu Y, Orre LM, Zhou Tran Y, et al. 2020. DEqMS: a method for accurate variance estimation in differential protein expression analysis. *Mol Cell Proteomics* 19(6):1047-1057.
- Ahlmann-Eltze C, Anders S. 2019. proDA: probabilistic dropout analysis for identifying differentially abundant proteins in label-free mass spectrometry. *bioRxiv* 661496 (preprint; cite `citation("proDA")`, never a journal).
- Choi M, Chang CY, Clough T, Broudy D, Killeen T, MacLean B, Vitek O. 2014. MSstats: an R package for statistical analysis of quantitative mass spectrometry-based proteomic experiments. *Bioinformatics* 30(17):2524-2526.
- Goeminne LJE, Gevaert K, Clement L. 2016. Peptide-level robust ridge regression improves estimation, sensitivity, and specificity in data-dependent quantitative label-free shotgun proteomics. *Mol Cell Proteomics* 15(2):657-668.
- Sticker A, Goeminne L, Martens L, Clement L. 2020. Robust summarization and inference in proteome-wide label-free quantification. *Mol Cell Proteomics* 19(7):1209-1219.
- Lazar C, Gatto L, Ferro M, Bruley C, Burger T. 2016. Accounting for the multiple natures of missing values in label-free quantitative proteomics data sets to compare imputation strategies. *J Proteome Res* 15(4):1116-1125.
- Stephens M. 2017. False discovery rates: a new deal. *Biostatistics* 18(2):275-294.
- Ebrahimpoor M, Goeman JJ. 2021. Inflated false discovery rate due to volcano plots: problem and solutions. *Brief Bioinform* 22(5):bbab053.

## Related Skills

- quantification - peptide-to-protein summarization, normalization, and IRS that produce the matrix this skill tests
- proteomics-qc - quality control and batch-effect assessment before testing
- protein-inference - razor/shared-peptide ambiguity that drives which protein group gets the quantity
- ptm-analysis - site-level differential testing for modified peptides
- differential-expression/de-results - analogous empirical-Bayes interpretation for RNA-seq DE
- data-visualization/volcano-and-ma-plots - volcano and MA plots of the result table
- pathway-analysis/go-enrichment - functional enrichment of the significant protein hit list
- machine-learning/biomarker-discovery - building predictive panels from differential proteins
- workflows/proteomics-pipeline - end-to-end pipeline that calls this skill as the testing stage
