---
name: bio-proteomics-differential-abundance
description: Tests for differentially abundant proteins between conditions with limma/DEqMS empirical-Bayes moderation, proDA/msqrob2/MSstats missingness modeling, and Python Welch+BH alternatives. Frames missing values as left-censored MNAR (model, do not impute), makes variance moderation the load-bearing step at n=3-5, and prefers feature/peptide-level testing. Use when identifying proteins with significant abundance changes between experimental groups. Summarization and normalization mechanics are proteomics/quantification; volcano and MA plots are data-visualization/volcano-and-ma-plots; pathway enrichment of the hit list is pathway-analysis/go-enrichment.
tool_type: mixed
primary_tool: limma
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: limma 3.62.2, DEqMS 1.24.0, proDA 1.20.0, msqrob2 1.14.1, QFeatures 1.16.0, MsCoreUtils 1.18.0, MSstats 4.14.2, ashr 2.2.63 (R 4.4.3 / Bioconductor 3.20), pandas 3.0.5, scipy 1.18.1, statsmodels 0.15.0 (checked 2026-09-15)

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

Scope: this skill owns the statistical TEST -- design/contrast construction, variance moderation, missingness handling, multiple-testing correction, minimum-fold-change testing, and fold-change shrinkage. Peptide-to-protein summarization and normalization mechanics -> proteomics/quantification. Volcano/MA plots -> data-visualization/volcano-and-ma-plots. Enrichment of the hit list -> pathway-analysis/go-enrichment. OUT OF SCOPE: how MaxLFQ/TMP/IRS produce the matrix (quantification); how to draw a volcano (data-visualization); single-sample (n=1) comparisons, and diagnosis or treatment decisions for an individual patient.

## The Single Most Important Modern Insight -- Model the Missingness, Moderate the Variance, Test at the Feature Level

1. **Missing values in label-free MS are left-censored MNAR -- missing BECAUSE the intensity is low -- and imputing them (especially Perseus/MaxQuant downshift) turns a detection limit into numbers the test treats as measurements. Downshift draws each missing value from a Gaussian with mean = mu - 1.8*sigma and SD = 0.3*sigma, where mu and sigma describe ALL proteins in that run (sigma is not the replicate SD). For an on/off protein (seen in all of group A, missing in all of B) the fold change is therefore set by the imputation constant and the protein's own intensity, not by biology -- the volcano-plot "anchor/wing" artifact (on/off proteins on a streak whose logFC tracks average intensity). The p-value carries no FDR guarantee: whether 0.3*sigma is narrower or wider than the real replicate SD decides whether imputation inflates false positives or only costs power. The honest statement is "undetected in group B", not "20x lower, p=1e-6". The correct approach is to MODEL the dropout in the likelihood -- proDA (probabilistic dropout), msqrob2, MSstats-AFT -- NOT fill it (Lazar 2016; Ahlmann-Eltze & Anders 2019).
2. **At the n=3-5 replicates proteomics actually uses, per-protein variance has only 2-4 residual df and is unusable raw -- variance moderation is the load-bearing element, not optional.** limma borrows a prior d0 across all proteins so a 4-replicate design tests on ~10 df instead of 6; `trend=TRUE` makes the prior a function of mean intensity (effectively mandatory for label-free, where a single global prior mis-calibrates FDR across the abundance range); `robust=TRUE` Winsorizes outlier variances (Phipson 2016). DEqMS makes the prior a function of PSM/peptide count and generally outperforms limma-trend when quantification depth varies across proteins (Zhu 2020).
3. **Feature/peptide-level modeling keeps information that summarize-then-test throws away -- but it is not free.** Summarizing first (one number per protein per run) discards the within-protein between-peptide variance and the correct degrees of freedom: 12 consistent peptides deserve a smaller SE than 12 disagreeing ones, but after summarization both look equally certain, and a protein with 30 observations looks as informative as one with 3. msqrob2/MSstats keep every peptide as a degree of freedom (Goeminne 2016; Sticker 2020; Choi 2014). Two calibrations from running all three on the same 4 v 4 peptide table: on clean, balanced data the feature-level fit and limma on the robustly summarized matrix returned the SAME 80 calls -- the gain is real when peptides within a protein disagree or coverage is unbalanced, not automatic. And the smaller SEs cut both ways: an uncorrected global offset that a protein-level test ignores becomes proteome-wide significance at the feature level, so a feature-level result is only readable after the centred-contrast check below.

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
| A peptide/precursor table is available, not just a protein matrix | msqrob2 (`QFeatures` + `msqrob`) or MSstats | the feature table is the only place the between-peptide spread still exists; both keep it |
| Outlier-peptide risk, unbalanced peptide coverage | msqrob2 `robustSummary` + `msqrob(robust = TRUE)` | Huber M-estimation downweights the bad peptide instead of letting it move the protein mean. On a clean, balanced 4 v 4 it buys nothing: msqrob2 and limma on the same summarized matrix returned the identical 80 calls. Escalate for the disagreement, not by default |
| Three or more groups, or group plus covariates, at the peptide level | msqrob2 `msqrobAggregate(..., ridge = TRUE)` | ridge shrinks effects estimated from few observations -- but it is REFUSED on a two-group design (needs >2 mean-model parameters); use `ridge = FALSE` there |
| Technical replicates, nested/repeated-measures, labeled (SRM/PRM/DIA) | MSstats (feature-level mixed model) | random effects capture run/subject structure summarize-then-test discards |
| Any feature-level route | check `median(log2FC) ~ 0` before reading the table | small feature-level SEs turn a normalization offset into proteome-wide significance; measured 21.2% realized FDR from a -0.18 offset |
| Batch present | batch as a covariate in the design (`~ batch + condition`) | `removeBatchEffect` is visualization-only; never feed its output to `lmFit` |
| Minimum biologically meaningful fold change | `treat()` + `topTreat()` (or SAM s0) | tests |log2FC|>c against the moderated null; a post-hoc FC+significance double filter inflates FDR |
| Large n (>10/group), Python-only | Welch t-test + BH | variance estimates reliable; no moderation needed at large n |

Default when uncertain: protein-level summary matrix at n=3-5 -> limma `eBayes(trend=TRUE, robust=TRUE)`; if PSM/peptide counts exist, escalate to DEqMS; if missingness is extensive and intensity-dependent, escalate to proDA.

## limma Workflow (R)

**Goal:** Identify differentially abundant proteins using moderated statistics that borrow information across all proteins.

**Approach:** Filter to proteins with enough valid values per group (rows with no or too few values give `NA` averages that stop `eBayes(trend = TRUE)`), build the design (batch as a covariate when present), fit the linear model, drop rows the design cannot estimate, then contrast, apply EB moderation with the intensity trend and robust fitting, and extract BH-corrected results. The per-condition count is NOT sufficient in a paired or blocked design (donor, subject or batch in the model): rows whose observed values fall in different blocks in the two conditions leave zero residual df or a partly `NA` coefficient vector, and limma tests them anyway on a contrast that is partly a block difference. Report proteins removed by either filter (e.g. undetected in one group) separately. Never feed `removeBatchEffect` output to `lmFit`.

```r
library(limma)

cond <- factor(sample_info$condition)
# Valid-value filter BEFORE lmFit: >= 2 values in every group (study choice; 3 of 4 is common)
n_valid <- sapply(levels(cond), function(g) rowSums(!is.na(protein_matrix[, cond == g, drop = FALSE])))
protein_matrix <- protein_matrix[apply(n_valid >= 2, 1, all), ]

design <- model.matrix(~0 + condition + batch, data = sample_info)  # batch in the model, not removed first
colnames(design)[seq_len(nlevels(cond))] <- levels(cond)

fit <- lmFit(protein_matrix, design)
# Estimability filter: the per-group count above does not make the contrast estimable under a blocked
# design ('Partial NA coefficients for N probe(s)'). Keep only fully estimated rows with residual df.
estimable <- fit$df.residual > 0 & rowSums(is.na(fit$coefficients)) == 0
fit <- fit[estimable, ]    # report the dropped rows; they are the non-estimable ones, not "not significant"

contrast_matrix <- makeContrasts(Treatment - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)  # trend mandatory for label-free; robust Winsorizes outliers

results <- topTable(fit2, coef = 1, number = Inf, adjust.method = 'BH')
# columns: logFC, AveExpr, t, P.Value, adj.P.Val, B  (adj.P.Val is the BH p; there is no $FDR)
```

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

## msqrob2 Workflow (R) -- Peptide Table In, Protein Calls Out

**Goal:** Test from the peptide/precursor table itself (MaxQuant `evidence.txt`, a DIA-NN report, a PSM table) so the between-peptide spread and the number of observations set the standard error, instead of collapsing to one number per protein per run first.

**Approach:** Build a `QFeatures` object from a wide peptide matrix, log-transform, keep peptides seen in at least 2 runs, aggregate to protein with `robustSummary` (Huber M-estimation, which is what downweights an outlier peptide), then fit `msqrob` -- ridge/robust regression with empirical-Bayes variance moderation -- and test the contrast. `hypothesisTest` writes its result into `rowData(pe[['protein']])`, one data frame per contrast. Two things the object does silently and you must undo: proteins with no observation in a condition come back with `adjPval = NA` rather than in a list, and any per-run normalization you apply here is applied to a *different peptide set in each run*, which is the offset trap below. Report the untestable and undetected proteins separately, as the limma section requires.

```r
library(QFeatures)
library(msqrob2)

# peptide_wide: one row per precursor; columns 'feature', 'protein', then one intensity column per run
runs <- sample_info$run
col_data <- data.frame(quantCols = runs, condition = factor(sample_info$condition),
                       sample = factor(runs), row.names = runs)  # quantCols column is required by readQFeatures
pe <- readQFeatures(assayData = peptide_wide, quantCols = runs, colData = col_data, name = 'peptideRaw')
pe <- zeroIsNA(pe, 'peptideRaw')
pe <- logTransform(pe, base = 2, i = 'peptideRaw', name = 'peptideLog')
rowData(pe[['peptideLog']])$nNonZero <- rowSums(!is.na(assay(pe[['peptideLog']])))
pe <- filterFeatures(pe, ~ nNonZero >= 2, keep = TRUE)  # keep=TRUE: the variable exists only on peptideLog

# undetected list, taken BEFORE aggregation: msqrob2 reports these as adjPval = NA, not as a list
cond <- colData(pe)$condition  # colData lives on the QFeatures object; pe[['peptideLog']]$condition is NULL
obs <- sapply(levels(cond), function(g)
  tapply(rowSums(!is.na(assay(pe[['peptideLog']])[, cond == g, drop = FALSE])),
         rowData(pe[['peptideLog']])$protein, sum))
undetected <- rownames(obs)[apply(obs, 1, min) == 0]  # report as "undetected in group X", never as a fold change

pe <- aggregateFeatures(pe, i = 'peptideLog', fcol = 'protein', name = 'protein',
                        fun = MsCoreUtils::robustSummary, na.rm = TRUE)
pe <- msqrob(pe, i = 'protein', formula = ~condition, robust = TRUE)
L <- makeContrast('conditionTreatment = 0', parameterNames = 'conditionTreatment')
pe <- hypothesisTest(pe, i = 'protein', contrast = L)

res <- rowData(pe[['protein']])$conditionTreatment  # columns: logFC, se, df, t, pval, adjPval (no adj.P.Val)
res$protein <- rownames(pe[['protein']])
untestable <- res$protein[is.na(res$adjPval)]       # report these; they are not "not significant"
res <- res[!is.na(res$adjPval), ]
```

To keep every peptide as its own degree of freedom instead of summarizing first, fit the mixed model over the peptide assay. `sample` must be a column of `colData` and `feature` a column of `rowData`:

```r
pe <- msqrobAggregate(pe, i = 'peptideLog', fcol = 'protein', name = 'proteinLmer',
                      formula = ~condition + (1 | sample) + (1 | feature), ridge = FALSE)
pe <- hypothesisTest(pe, i = 'proteinLmer', contrast = L)
```

`ridge = TRUE` needs MORE than two parameters in the mean model, so it is refused outright on a plain two-group comparison ("The mean model must have more than two parameters for ridge regression"); use `ridge = FALSE`, or drop the intercept with `~ -1 + condition` as the error message itself suggests. Ridge is for three or more groups, or a group factor plus covariates.

## MSstats Workflow (R) -- Feature-Level Mixed Model

**Goal:** Test a feature-level design with run/subject structure -- technical replicates, nested or repeated measures, SRM/PRM/DIA -- from the search engine's own output tables.

**Approach:** Convert with the vendor-specific importer, summarize with `dataProcess`, then test the contrast with `groupComparison`. The contrast matrix is columns-by-condition and its `dimnames` must match `levels(proc$ProteinLevelData$GROUP)` exactly. Keep `MBimpute = FALSE` unless you want the AFT censored imputation: it is the only place MSstats imputes, and the Skill's position is to model dropout rather than fill it. On a 4 v 4 label-free set the two settings differed by 0.4 percentage points of realized FDR, so the AFT imputation is not buying accuracy either. Proteins present in only one condition come back with an infinite `log2FC` and `issue == 'oneConditionMissing'`; split them out and report them as undetected.

```r
library(MSstats)

input <- MaxQtoMSstatsFormat(evidence = evidence, proteinGroups = protein_groups,
                             annotation = annotation, use_log_file = FALSE)  # annotation: Raw.file, Condition, BioReplicate, IsotopeLabelType
proc <- dataProcess(input, normalization = 'equalizeMedians', summaryMethod = 'TMP',
                    censoredInt = 'NA', MBimpute = FALSE, use_log_file = FALSE)

contrast <- matrix(c(-1, 1), nrow = 1,
                   dimnames = list('Treatment-Control', c('Control', 'Treatment')))  # order = levels(GROUP)
res <- groupComparison(contrast.matrix = contrast, data = proc, use_log_file = FALSE)$ComparisonResult
res$Protein <- as.character(res$Protein)

undetected <- res[res$issue %in% 'oneConditionMissing', ]   # infinite log2FC; report as undetected, not as a ratio
tested <- res[is.finite(res$log2FC) & !is.na(res$adj.pvalue), ]
# columns: Protein, Label, log2FC, SE, Tvalue, DF, pvalue, adj.pvalue, issue (adj.pvalue is the BH p)
```

### Check the contrast is centred before reading any feature-level result

**Goal:** Catch the one failure that turns a feature-level test anticonservative across the whole table.

**Approach:** A feature-level route has far smaller standard errors than a protein-summary test -- that is the point of it -- so a global between-condition offset of 0.1-0.2 log2 that limma would never call becomes significant for hundreds of proteins at once. The offset comes from per-run median normalization applied to the peptide table, because each run's median is taken over a *different* set of detected peptides. Measure it: the median log2FC over all tested proteins must be ~0 unless you genuinely expect most of the proteome to move. If it is not, the normalization is the suspect, not biology -- fix it upstream in proteomics/quantification (normalize on features present in every run, or at the protein level after summarization) rather than by re-centring the p-values here.

```r
offset <- median(tested$log2FC, na.rm = TRUE)          # msqrob2: median(res$logFC)
if (abs(offset) > 0.05) stop(sprintf(
  'median log2FC = %+.3f: the contrast is not centred. Re-normalize (proteomics/quantification) before reading this table.',
  offset))
```

## Python Workflow

**Goal:** Run the full pipeline in Python when no R is available and n is large enough that moderation is unnecessary.

**Approach:** Log2-transform, median-normalize, run per-protein Welch t-tests, apply Benjamini-Hochberg. Return the untestable proteins (fewer than 2 values in a group) alongside the results, as the limma section requires -- they are not "not significant". This has NO variance moderation and should not be used at n=3-5 -- escalate to limma/DEqMS for small n.

```python
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

def preprocess(intensities):
    log2_data = np.log2(intensities.replace(0, np.nan))  # zeros -> NaN to avoid -inf
    sample_medians = log2_data.median(axis=0)
    return log2_data - sample_medians + sample_medians.median()

def differential_abundance(normalized, case_cols, ctrl_cols):
    rows, untestable = [], []
    for protein in normalized.index:
        case, ctrl = normalized.loc[protein, case_cols].dropna(), normalized.loc[protein, ctrl_cols].dropna()
        if len(case) >= 2 and len(ctrl) >= 2:
            _, pval = stats.ttest_ind(case, ctrl, equal_var=False)  # Welch; scipy defaults to Student's True
            rows.append({'protein': protein, 'log2fc': case.mean() - ctrl.mean(), 'pvalue': pval})
        else:  # never drop these silently: proteins undetected in one group are often the largest real changes
            untestable.append({'protein': protein, 'n_case': len(case), 'n_ctrl': len(ctrl)})
    if not rows:
        raise ValueError('No protein has >= 2 non-missing values in both groups; a two-sample test is not possible')
    df = pd.DataFrame(rows)
    df['padj'] = multipletests(df['pvalue'], method='fdr_bh')[1]  # default is Holm-Sidak; pass fdr_bh explicitly
    return df, pd.DataFrame(untestable, columns=['protein', 'n_case', 'n_ctrl'])  # report the second table too
```

## Fold-Change Reporting

**Goal:** Hand the right effect estimate to the right consumer.

**Approach:** Report the RAW fold change (the best unbiased point estimate) for GSEA/pathway ranking and meta-analysis -- those need the full continuous distribution or FC+SE pairs. Apply shrinkage (ashr) only when recovering "which proteins truly changed and by how much"; it fits a mixture prior with a point mass at zero and shrinks uncertain effects smoothly toward zero. This is preferred over hard-thresholding (zeroing FCs at padj 0.05), which creates an arbitrary step function. No mature Python ashr equivalent exists.

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
**Mechanism:** Imputed values come from one Gaussian set by the run-wide intensity distribution -> the fold change of an on/off protein is fixed by the downshift constant, and the imputed within-group SD (0.3 x run-wide SD) is unrelated to the real replicate SD, so the test is anticonservative or merely underpowered depending on the data.
**Symptom:** Volcano "anchor/wing" -- on/off proteins on a streak whose logFC tracks average intensity; p-values and FDR without a guarantee.
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

### Per-run median normalization under a feature-level test
**Trigger:** `dataProcess(normalization = 'equalizeMedians')`, `normalize(method = 'center.median')` or any per-run median centring applied to the PEPTIDE table, then msqrob2/MSstats.
**Mechanism:** Each run's median is computed over the peptides detected in THAT run, and the detected set differs run to run. Equalizing those medians therefore transfers a detection-composition difference into a uniform between-condition offset in every protein. A protein-summary test would mostly shrug it off; a feature-level test has standard errors small enough to call it significant, so the offset becomes a proteome-wide false-positive band.
**Symptom:** `median(log2FC)` over all tested proteins is far from 0 while the raw peptide table's is ~0; the extra calls are low-|logFC|, low-SE proteins; the true hits are recovered either way.
**Fix:** Check `median(log2FC)` before reading the table. Normalize on features present in every run, or after summarization at the protein level (proteomics/quantification owns the mechanics). Do not "fix" it by re-centring p-values.

### FC + significance double filter
**Trigger:** `abs(logFC) > 1 & adj.P.Val < 0.05` applied after the test.
**Mechanism:** |logFC| is large for a true effect OR a large SE; filtering on both the FC and the p (both depend on SE) selects high-variance nulls (collider effect).
**Symptom:** Realized FDR above nominal with no guarantee -- above 50% for top-ranked lists with many small effects (Ebrahimpoor & Goeman 2021), much lower when true effects are large.
**Fix:** `treat()`+`topTreat()` or SAM s0, which sit inside the statistic before selection.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| n=3-5 replicates -> 2-4 residual df | -- | raw per-protein variance unusable; moderation is mandatory, not optional |
| limma adds prior d0 (~4) df | Ritchie 2015 | a 4-replicate design tests on ~10 df vs 6; the borrowed df is the benefit |
| downshift mean = mu - 1.8*sigma, SD = 0.3*sigma | Perseus default | sigma is the run-wide SD across proteins, not the replicate SD; the imputed FC is set by the constant, and 0.3*sigma can be wider or narrower than real replicate spread |
| `trend=TRUE` effectively mandatory for label-free | Ritchie 2015 | a single global prior mis-calibrates FDR across abundance |
| min-FC floor log2(1.2) (1.2-fold) via treat() | -- | example floor; common alternatives 1.5-fold (~0.58) or 2-fold (1.0); set by biology, tested against the moderated null |
| BH adjusted p < 0.05 | Benjamini-Hochberg | controls FDR over the WHOLE rejection set, not subsets carved out afterward |
| DEqMS multi-batch TMT: minimum count across batches | Zhu 2020 | the bottleneck batch sets the realized precision |
| realized FDR > 50% from FC+significance double filter | Ebrahimpoor & Goeman 2021 | top-100 at n=12 exceeded 50% FDR at nominal 5%; regime-specific (many small effects), not a general rate |
| feature-level route: `abs(median(log2FC))` <= 0.05 | measured | on one 4 v 4 peptide set the realized FDR tracked the offset monotonically: -0.005 -> 0.0%, -0.107 -> 7.0%, -0.184 -> 21.2%, -0.202 -> 28.6%, all at nominal 5%; the true-hit count did not change |
| msqrob2 ridge needs > 2 mean-model parameters | msqrob2 1.14.1 | a two-group `~condition` mean model has 2, so `ridge = TRUE` errors; `ridge = FALSE` or `~ -1 + condition` |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `results$FDR` is NULL | limma `topTable`/`topTreat` have no `$FDR` column | use `adj.P.Val` (BH-adjusted p) |
| `topTreat` row has no `B` | `topTreat` omits `B` (a `topTable` column) | read `logFC, AveExpr, t, P.Value, adj.P.Val` |
| FDR mis-calibrated across abundance | `eBayes` with `trend=FALSE` on intensity data | `eBayes(fit, trend = TRUE, robust = TRUE)` |
| min-FC test inflates FDR | `topTable(lfc=...)` or post-hoc volcano double filter | `treat(fit, lfc=log2(1.2), trend=TRUE, robust=TRUE)` then `topTreat()` |
| treat() list moderated without the intensity trend | `treat()` re-estimates the prior with `trend=FALSE, robust=FALSE` by default | pass `trend = TRUE, robust = TRUE` to `treat()` |
| `eBayes`: `prior.weights contain NA values` | rows with no valid value (MaxQuant all-zero LFQ rows) or too few per group | valid-value filter before `lmFit` |
| `lmFit`: `Partial NA coefficients for N probe(s)`; `stopifnot(all(fit2$df.residual > 0))` fires | a paired/blocked design (donor, subject, batch) where a row's values sit in different blocks per condition; the per-condition valid-value filter does not catch this | `fit <- fit[fit$df.residual > 0 & rowSums(is.na(fit$coefficients)) == 0, ]` after `lmFit` |
| DEqMS warning `longer object length is not a multiple of shorter object length` | rows with NA sigma / zero residual df reach `spectraCounteBayes` | filter before `lmFit`; require `fit2$df.residual > 0` |
| proDA: `object 'conditionControl' not found` | intercept design with `reference_level`: coefficients are `Intercept`, `conditionTreatment` | `test_diff(fit, 'conditionTreatment')` |
| `makeContrasts`: `object 'DrugB' not found` | design rename hard-coded to two groups | `colnames(design)[seq_len(nlevels(cond))] <- levels(cond)` |
| Python `ValueError: No protein has >= 2 non-missing values` | a group has one sample, or no protein is observed twice per group | not a two-sample test; do not run DA on n=1 |
| anticonservative p after batch correction | `removeBatchEffect` output fed to `lmFit` | put batch in the design: `~ batch + condition` |
| DEqMS columns missing | forgot `fit$count` or read limma columns | set `fit$count`, run `spectraCounteBayes`, read `sca.adj.pval` from `outputResult` |
| Student's t instead of Welch | `scipy.stats.ttest_ind` defaults `equal_var=True` | pass `equal_var=False` |
| p-values look like Holm-Sidak | `statsmodels` `multipletests` defaults to `'hs'` | pass `method='fdr_bh'` |
| volcano "anchor/wing" streaks | downshift/imputation feeding the test | model dropout (proDA/msqrob2/MSstats-AFT); report on/off proteins as undetected |
| hundreds of low-|logFC| calls from msqrob2/MSstats and `median(log2FC)` far from 0 | per-run median normalization on the peptide table offsets the whole contrast | re-normalize on complete-case features or after summarization; see the centred-contrast check |
| `readQFeatures`: `'colData' must contain a column called 'quantCols'` | QFeatures 1.16 requires the run names in `colData$quantCols` | `data.frame(quantCols = runs, condition = ..., row.names = runs)` |
| msqrob2: `The mean model must have more than two parameters for ridge regression` | `ridge = TRUE` on a two-group design | `ridge = FALSE`, or drop the intercept (`~ -1 + condition`) |
| msqrob2: `Variable sample is not found in coldata or rowdata` | `msqrobAggregate` formula names a term that is in neither | add `sample` to `colData`; `feature` must be a `rowData` column of the peptide assay |
| `pe[['peptideLog']]$condition` is `NULL` | in QFeatures 1.16 `colData` is held on the QFeatures object, not on each assay | `colData(pe)$condition` (or `getWithColData(pe, i)` for a standalone assay) |
| msqrob2 results have no `adj.P.Val`; some rows are all `NA` | `hypothesisTest` writes `logFC, se, df, t, pval, adjPval` into `rowData`, and returns `NA` for proteins it cannot fit | read `adjPval`; report the `is.na(adjPval)` rows separately, not as non-significant |
| MSstats rows with infinite `log2FC` | `issue == 'oneConditionMissing'` | split them out and report as undetected; never as a ratio |
| MaxQuant `evidence.txt`: the `Reverse`/`Potential contaminant` filter drops every row | empty flag columns are read as logical `NA`, and `NA != '+'` is `NA` | filter with `!(ev$Reverse %in% '+')`, which is `FALSE` for `NA` |

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
