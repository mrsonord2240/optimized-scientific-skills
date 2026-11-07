# Background

Not needed to run an analysis. This explains choices the routes and scripts already make. Scope of the Skill: the statistical test on a log2 normalized matrix or a peptide table. Peptide-to-protein summarization and normalization, volcano and MA plots, and enrichment of the hit list belong to other Skills (quantification, data-visualization, pathway-analysis). Single-sample comparisons and decisions about an individual patient are out of scope.

## Model the missingness, moderate the variance, test at the feature level

1. Missing values in label-free MS are left-censored and not at random: missing because the intensity is low. Perseus/MaxQuant downshift draws each missing value from a Gaussian with mean mu - 1.8 sigma and SD 0.3 sigma, where mu and sigma describe all proteins in that run. For an on/off protein (seen in all of group A, missing in all of B) the fold change is then set by the imputation constant, giving the volcano "anchor/wing" streak. Whether the 0.3 sigma SD is narrower or wider than the real replicate SD decides whether imputation inflates false positives or only costs power. The honest statement is "undetected in group B". proDA, msqrob2 and MSstats-AFT model the dropout instead.
2. At 3 to 5 replicates a per-protein variance has 2 to 4 residual df. limma borrows a prior d0 across proteins, so a 4-replicate design tests on about 10 df instead of 6. `trend = TRUE` ties the prior to mean intensity (needed for label-free); `robust = TRUE` Winsorizes outlier variances. DEqMS ties the prior to PSM or peptide count and does better than limma-trend when quantification depth varies.
3. Summarizing first discards the between-peptide variance and the true degrees of freedom: 12 consistent peptides and 12 disagreeing ones look equally certain afterwards. msqrob2 and MSstats keep every peptide. On one clean balanced 4 v 4 peptide table msqrob2 and limma on the robustly summarized matrix made the same 80 calls; the gain appears when peptides disagree or coverage is unbalanced.

## Tools

| Tool | Mechanism | Use |
|------|-----------|-----|
| limma | EB moderated t; prior d0 blended with the per-protein variance | protein-level summaries, small n |
| DEqMS | prior variance is a loess of log-variance on log2(count) | TMT (PSM count) and label-free DDA (peptide count) |
| proDA | probabilistic dropout: left-censored values integrated under a per-sample sigmoid; EB on location and variance | label-free DDA with many intensity-dependent missing values |
| msqrob2 | peptide-level robust ridge on a `QFeatures` object: Huber M-estimation downweights outlier peptides, ridge shrinks effects from few observations (needs more than 2 mean-model parameters), EB variance | outlier-peptide or unbalanced-coverage risk |
| MSstats | feature-level linear mixed model (group fixed; feature and run/subject random) | SRM/PRM/DIA, technical replicates, nested or repeated measures |
| Welch + BH | per-protein two-sample t, BH | more than 10 per group, Python only |
| ashr | mixture prior with a point mass at zero; posterior means shrink uncertain effects | which proteins truly changed and by how much |

## Measurements behind the scripts

- On a synthetic 6 v 6 set with a 0.34 log2 loading offset, limma on the un-normalized matrix made 135 calls at 30% false discoveries; after median centring, 109 calls at 7%. The `limma_de.R` median log2FC stop is 0.1.
- proDA on 4 v 4: 0 of 11 on/off proteins called, best adjusted p 0.17. It is honest but underpowered for them.
- Per-run median normalization of the peptide table gave 20.8% to 27.9% realized FDR at nominal 5% (median log2FC -0.18 to -0.20, residual SD 0.76 of its un-normalized value). MSstats `equalizeMedians`: median log2FC -0.184, median null SE 0.119, 101 calls with 21 false positives; `normalization = FALSE`: +0.008, SE 0.204, 79 calls, 0 false positives. On MSstats, MBimpute TRUE vs FALSE differed by 0.4 percentage points of realized FDR.
- Decomposition on one 4 v 4 peptide set, msqrob2:

| variant | median log2FC | null-protein median SE | calls | false positives |
|---|---|---|---|---|
| no normalization | -0.005 | 0.206 | 80 | 0 |
| uniform -0.20 injected into every treatment run | -0.205 | 0.206 | 79 | 0 |
| each run centred on its own condition's median | -0.016 | 0.116 | 81 | 1 |
| per-run median over all detected peptides | -0.199 | 0.116 | 111 | 31 (27.9%) |

- `OFFSET_MAX` 0.05 and `SD_RATIO_MIN` 0.85 in the feature-level checks are trip-wires read off that one set, not distributional bounds. Realized FDR was 0.0% at offset -0.005, 8.0% at -0.107, 27.9% at -0.199. The offset stop would block a harmless uniform shift; the within-condition variant passes the offset check and is caught only by the SD ratio.

## Failure mechanisms

- Imputation feeding a variance-based test: the on/off fold change is fixed by the downshift constant and the imputed SD is unrelated to the replicate SD. kNN or mean imputation is mean-reverting, so truly low values are pulled up and down-regulation is compressed; it is valid only under random missingness.
- `removeBatchEffect` before `lmFit`: subtracts the batch component without propagating uncertainty, so residual variance is understated and EB df inflated; if batch is confounded with biology it deletes real signal.
- `eBayes(trend = FALSE)` on intensities: one global prior over-shrinks high-abundance and under-shrinks low-abundance proteins.
- Wrong DEqMS count (razor+unique vs MS2 PSMs; total vs minimum across batches): the prior is fit on the wrong precision proxy.
- proDA on random missingness: the censored-dropout model is mis-specified.
- FC + p double filter: BH controls the rate for FC = 0, not for the FC threshold; conditioning on both selects high-variance nulls. Inflation is above 50% for top-ranked lists with many small effects, far smaller when true effects are large. `treat()` tests against the moderated null instead and re-estimates its prior, so `trend` and `robust` are passed again (they default to FALSE).
- Shrinkage: raw fold change is the best unbiased point estimate, so GSEA ranking and meta-analysis use it (pool raw FCs with their SEs and shrink afterwards). Shrinkage replaces hard-thresholding, which creates an arbitrary step.
- msqrob2 and DEqMS: DEqMS columns are missing without `fit$count`; msqrob2 `ridge = TRUE` needs more than 2 mean-model parameters.

## Thresholds

| Threshold | Note |
|-----------|------|
| `fold` floor | 1.2-fold is an example; 1.5-fold (log2 0.58) and 2-fold (1.0) are common. Set by biology. |
| BH adjusted p < 0.05 | Controls FDR over the whole rejection set, not subsets carved out afterwards. |
