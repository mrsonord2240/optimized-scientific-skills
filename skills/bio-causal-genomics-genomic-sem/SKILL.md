---
name: bio-causal-genomics-genomic-sem
description: Fits structural equation models to GWAS summary statistics using GenomicSEM (Grotzinger 2019), including common-factor models, confirmatory factor models, ESEM, common-factor GWAS with Q_SNP heterogeneity, multivariate Wald tests, and stratified GenomicSEM partitioned heritability. Reconciles results against MTAG multi-trait analysis. Handles sample overlap via the LDSC sampling-covariance matrix, identifies and resolves Heywood cases, and verifies model fit with CFI / RMSEA. Use when modeling latent genetic architecture across correlated traits, running multivariate GWAS on a shared factor, distinguishing factor-mediated from trait-specific SNP effects, or comparing GenomicSEM common-factor results against MTAG when both depend on accurate sampling covariance.
tool_type: r
primary_tool: GenomicSEM
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: GenomicSEM 0.0.5 (GitHub `GenomicSEM/GenomicSEM`) + **lavaan 0.6.19 pinned** (see Tool Installation for the install command and why), LDSC (Python 3; `CBIIT/ldsc` commit `1f09cf0` -- see Tool Installation), baselineLD_v2.2 annotations (alkesgroup.broadinstitute.org/LDSCORE), MTAG 1.0.8+ (Python; `JonJala/mtag`), R 4.4+.

**lavaan is capped, not floored: use 0.6.19.** On lavaan >= 0.7.0, GenomicSEM 0.0.5's internal reorder-step `sem()` calls omit the `ordered = FALSE` that lavaan now requires for DWLS on continuous data, so `usermodel()`, `commonfactorGWAS()` and `userGWAS()` crash under both estimators, and no `estimation=` choice avoids it (symptoms in Common Errors). Upgrading GenomicSEM does not help: the calls are unchanged at GitHub HEAD `6b65ca5` (2026-08-26). GenomicSEM 0.0.5 + lavaan 0.6.19 runs all four core functions under DWLS and ML, recovering planted loadings and factor correlations on synthetic inputs.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('GenomicSEM')`, `packageVersion('lavaan')` then `?ldsc`, `?commonfactor`, `?usermodel`, `?commonfactorGWAS`
- Python (LDSC, MTAG): `<tool>.py -h` and inspect the source under `ldsc/` or `mtag/`

GenomicSEM is GitHub-only (never on CRAN). If `ldsc()` or `usermodel()` throws an error about lavaan syntax or non-positive-definite covariance -- and it is not the lavaan-version crash above -- introspect the installed API (`getMethod('ldsc')`) and adapt rather than retrying.

## Scope Boundary

GenomicSEM models latent genetic architecture across GWAS summary statistics at the population level; it does not predict or diagnose outcomes for any individual. Decline and redirect individual-level polygenic-score / diagnostic requests (e.g. "does this patient's PGS mean they will develop the disorder?") to a qualified clinician or genetic counselor -- this Skill's output is never an individual risk estimate.

# Genomic SEM

**"Model the latent genetic architecture across several correlated GWAS"** -> Treat each GWAS as a measured indicator of one or more latent genetic factors and fit a structural equation model to the LDSC-derived genetic covariance matrix S and its sampling covariance V (Grotzinger 2019 Nat Hum Behav 3:513). The framework extends naturally to a multivariate GWAS in which a SNP is regressed on a latent factor (common-factor GWAS), with Q_SNP testing whether the SNP effect is homogeneous across factor loadings. Sample overlap between input GWAS is absorbed by the off-diagonals of V; ignoring V inflates Type-I.

- R: `GenomicSEM::ldsc()` produces the (S, V) covariance pair from munged sumstats
- R: `GenomicSEM::commonfactor()` fits a single-factor CFA across all traits in S
- R: `GenomicSEM::usermodel()` fits an arbitrary lavaan-syntax model
- R: `GenomicSEM::commonfactorGWAS()` runs SNP -> factor multivariate GWAS with Q_SNP
- R: `GenomicSEM::userGWAS()` runs arbitrary multivariate SNP regression with per-path Wald tests and a model chi-square (no `Q_pval` column)
- Python (alternative): `mtag.py --sumstats t1,t2,t3 --out mtag_out` (multi-trait power boost on individual traits)

## Statistical Model Taxonomy

| Method | Latent structure | Min traits | SNP-level test | Strength | Fails when |
|--------|------------------|-----------|----------------|----------|------------|
| Common-factor CFA (Grotzinger 2019) | Single F loading all traits | 3 | None (model-fit only) | Tests whether shared variance is unidimensional | Heterogeneous architecture; CFI < 0.9; near-zero loadings |
| User-specified CFA (`usermodel`) | Pre-specified lavaan syntax | 3 | None | Confirmatory; arbitrary structure | Misspecified model; identification under-determined |
| ESEM | Exploratory rotation; cross-loadings allowed | 6+ | None | When factor count and structure unknown | Few traits; collinear traits; rotation arbitrary |
| Common-factor GWAS (`commonfactorGWAS`) | SNP -> F -> trait1..k | 3 | Wald on F + Q_SNP heterogeneity | Discovers SNPs acting via the common factor; flags Q_SNP outliers | Q_SNP-significant SNPs not interpretable as factor SNPs |
| User GWAS (`userGWAS`) | Arbitrary SNP-path lavaan | 3 | Wald per path + model chi-square | Tests SNP on any specified path | Highly parameterized models lose power |
| Multivariate Wald test | Joint test across SNP -> trait paths | 2+ | Joint chi-square | Boost power when SNP affects multiple traits | Heterogeneous SNP effects collapse joint test |
| Stratified GenomicSEM (Grotzinger AD et al 2022 Nat Genet 54:548) | Factor model with sLDSC-partitioned annotations | 3 | Per-annotation factor tau | Localizes heritability of the factor to functional categories | Same sLDSC failure modes (small annotation, collinearity) |
| MTAG (Turley 2018 Nat Genet 50:229) | Empirical-Bayes shrinkage across correlated traits | 2 | Per-trait shrunk z-score | Boosts marginal power for any input trait | MaxFDR > 5% indicates heterogeneity violates MTAG assumption |

Methodology evolves; verify the current Grotzinger 2023+ tutorials at `github.com/GenomicSEM/GenomicSEM/wiki` before locking a method. ESEM rotation choice (geomin vs target rotation) is an active area; report sensitivity to rotation.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Multi-trait GWAS power boost for one focal trait | MTAG | Optimized for per-trait marginal power |
| Common-factor architecture hypothesized | GenomicSEM `commonfactorGWAS` | Tests SNP -> factor; reports Q_SNP heterogeneity |
| Heterogeneous architecture (>1 latent factor) | ESEM, then confirmatory `usermodel` (`references/advanced-models.md`) | Exploratory first, then confirm |
| Confirming a pre-specified factor structure | `usermodel` with lavaan syntax | Confirmatory factor analysis |
| Partition heritability of factor across annotations | Stratified GenomicSEM (`references/stratified-genomicsem.md`) | Combines sLDSC + factor model |
| Mediation in a SEM framework | `usermodel` with indirect path | Path coefficients + delta-method SE |
| Sample overlap unknown or any-overlap suspected | Always use `ldsc()` output as input | V matrix off-diagonals absorb overlap |
| Cross-ancestry common-factor analysis | Run per-ancestry, compare loadings; no published cross-ancestry SEM as of 2026 | Method not yet validated for mixed-ancestry V |
| Single biobank for all traits (e.g., UKB only) | GenomicSEM with `ldsc()`; the V matrix will reflect overlap | Equivalent to one-sample MR -- the V matrix is the correction |
| Comparing GenomicSEM and MTAG on the same traits (`references/mtag-comparison.md`) | Run both; compare top hits + heterogeneity | Concordance increases confidence; divergence flags heterogeneity |

## Per-Method Failure Modes

Heywood case, forced common factor on heterogeneous traits and non-positive-definite V: `references/failure-modes.md`. MTAG MaxFDR > 5%: `references/mtag-comparison.md`. Sample overlap and Q_SNP below.

### Sample overlap mis-specified

**Trigger:** Using LDSC intercept manually or supplying covariance from non-`ldsc()` source.

**Mechanism:** GenomicSEM's `ldsc()` function returns a list with `S` (genetic covariance) AND `V` (sampling covariance of the lower-triangle of S). The V off-diagonals capture sample overlap via cross-trait LDSC intercept. Skipping V and supplying only S treats all inputs as independent samples; Type-I error inflates because the sampling distribution under H0 is wrong.

**Symptom:** SE on factor loadings far too small; many SNPs significant in common-factor GWAS that don't replicate; comparison to MTAG shows disagreement consistent with overlap.

**Fix:** Always pass the full output of `ldsc()` -- both S and V -- to `commonfactor()`, `usermodel()`, and `commonfactorGWAS()`. Never construct S manually from rg estimates.

### Q_SNP not reported in commonfactorGWAS

**Trigger:** Running `commonfactorGWAS()` and reporting only the factor p-value per SNP.

**Mechanism:** Q_SNP tests heterogeneity of the SNP's effect across factor loadings (Grotzinger 2019 supplement). A SNP with significant Q_SNP violates the common-factor assumption: its effect is NOT mediated by the factor, and the factor estimate is meaningless for that SNP.

**Symptom:** Top "common-factor SNPs" are dominated by trait-specific effects; replication in independent cohorts is poor for SNPs with high Q_SNP.

**Fix:** Always report Q_SNP p-value alongside the factor p-value. Flag SNPs with Q_SNP p < 0.05 / N_factor_SNPs (Bonferroni for the discovered set; see Quantitative Thresholds) as architecture-violating and exclude from "common-factor SNP" claims. Re-fit those SNPs in `userGWAS()` with separate paths to each trait.

## Model Fit Diagnostics

| Index | Acceptable | Good | Source |
|-------|-----------|------|--------|
| CFI | >= 0.90 | >= 0.95 | Hu & Bentler 1999 Struct Equ Model 6:1 |
| TLI / NNFI | >= 0.90 | >= 0.95 | Hu & Bentler 1999 |
| RMSEA | <= 0.08 | <= 0.05 | Hu & Bentler 1999 |
| SRMR | <= 0.08 | <= 0.05 | Hu & Bentler 1999 |
| chi-square p-value | (less informative at large N) | n/a | Penalize for N inflation |

AIC / BIC are used for nested-model comparison (lower is better); only compare nested models fit on the same S.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Fit indices (CFI, RMSEA, SRMR) | see Model Fit Diagnostics | Report RMSEA 0.05 - 0.08 as "adequate", not "good" |
| Q_SNP p < 0.05 / N_SNP_factor | Grotzinger 2019 Nat Hum Behav 3:513 | Bonferroni for heterogeneity among factor-significant SNPs |
| MTAG MaxFDR < 5% | Turley 2018 Nat Genet 50:229 | Above this, MTAG marginal trait results invalid |
| Per-trait LDSC mean chi-square > 1.02 | LDSC documentation | Below this, V entries too noisy; factor SE inflated |
| Standardized loading 0.3 - 0.9 typical | SEM conventions | < 0.3 trait loads weakly; > 0.95 may indicate over-fit / collinearity |
| Min 3 traits for common factor | SEM identification | Single factor with k traits has k(k+1)/2 moments; needs k>=3 to identify |

## Standard Workflow

**Goal:** Fit a common-factor model across correlated GWAS and run a multivariate GWAS on the factor with Q_SNP.

**Approach:** Munge sumstats -> LDSC for (S, V) -> common-factor CFA -> inspect fit -> prepare SNPs -> common-factor GWAS -> report factor effects with Q_SNP flags.

```r
library(GenomicSEM)

# Step 1: Munge sumstats (one-time; produces .sumstats.gz files). GenomicSEM::munge()
# or LDSC's own munge_sumstats.py both work; either way, verify HapMap3-alignment first.
files <- c('raw/trait1.txt', 'raw/trait2.txt', 'raw/trait3.txt')
hm3 <- 'w_hm3.snplist'  # HapMap3 SNP list
trait_names <- c('trait1', 'trait2', 'trait3')
N <- c(150000, 200000, 175000)
munge(files = files, hm3 = hm3, trait.names = trait_names, N = N)

# Step 2: LDSC produces both S (genetic covariance) and V (sampling covariance)
traits <- c('trait1.sumstats.gz', 'trait2.sumstats.gz', 'trait3.sumstats.gz')
ldsc_results <- ldsc(
    traits = traits,
    sample.prev = c(0.5, 0.5, NA),    # case prevalence; NA for continuous
    population.prev = c(0.05, 0.05, NA),
    ld = 'eur_w_ld_chr/',
    wld = 'eur_w_ld_chr/',
    trait.names = trait_names
)
# ldsc_results$S = genetic covariance; ldsc_results$V = sampling covariance

# Step 3a: Common-factor CFA via DWLS
cf_fit <- commonfactor(covstruc = ldsc_results, estimation = 'DWLS')
print(cf_fit$modelfit)  # CFI, RMSEA, SRMR, chi-square
print(cf_fit$results)   # loadings + SEs

# Step 3b: Alternative -- user-specified two-factor model (>=3 indicators per factor)
# Identification rule: each factor needs >= 3 indicators OR one anchor loading fixed
# to 1 plus factor variance free. A factor with a single indicator is NOT identified.
model_syntax <- '
    F1 =~ NA*trait1 + trait2 + trait3
    F2 =~ NA*trait4 + trait5 + trait6
    F1 ~~ 1*F1
    F2 ~~ 1*F2
    F1 ~~ F2
'
user_fit <- usermodel(covstruc = ldsc_results, model = model_syntax, estimation = 'DWLS')
```

`estimation = 'DWLS'` (diagonally weighted least squares) is the default and is required when V is large; `'ML'` is faster but assumes a known V and can produce wrong SE under sample overlap. For the one exception, Q_SNP classification, see Common-Factor GWAS with Q_SNP.

**Result column names differ by function -- extract by name only after `names()`** (checked on GenomicSEM 0.0.5 + lavaan 0.6.19):

| Function | `$results` columns |
|----------|--------------------|
| `commonfactor()` | `lhs, op, rhs, Unstandardized_Estimate, Unstandardized_SE, Standardized_Est, Standardized_SE, p_value` |
| `usermodel()` | `lhs, op, rhs, Unstand_Est, Unstand_SE, STD_Genotype, STD_Genotype_SE, STD_All, p_value` |

The standardized loading is `Standardized_Est` in one and `STD_Genotype` in the other; `cf_fit$results[, 'STD_Genotype']` on a `commonfactor()` fit fails. `commonfactorGWAS()` and `userGWAS()` return per-SNP data frames with no standardized column (columns under their sections below).

## Reference Files

Read the file only when the task needs it; the rest of this skill is enough for a common-factor GWAS.

| Need | File |
|------|------|
| ESEM (exploratory factor structure), `userGWAS()` custom SNP path models, higher-order / bifactor / p-factor models | `references/advanced-models.md` |
| Run MTAG beside GenomicSEM, compare it with `commonfactorGWAS` (property table), reconcile the two when they disagree, or MTAG MaxFDR > 5% | `references/mtag-comparison.md` |
| Heywood case (negative residual variance, loading > 1), poor common-factor fit (CFI < 0.9), or non-positive-definite V | `references/failure-modes.md` |
| Runtimes and cluster settings for genome-wide runs, or installing the LDSC / MTAG Python tools | `references/runtime-and-python-tools.md` |
| Stratified GenomicSEM (`enrich()`, partitioned heritability of the factor) | `references/stratified-genomicsem.md` |

## Common-Factor GWAS with Q_SNP

**Goal:** Identify SNPs that affect the latent factor and flag SNPs whose effect is heterogeneous across loadings.

**Approach:** Build SNP-by-trait effect-and-SE matrix via `sumstats()`, then fit the SNP-augmented model genome-wide via `commonfactorGWAS()`. Report both factor p-value and Q_SNP p-value per SNP.

`sumstats()` flips effect signs based on the reference panel A1/A2 to enforce consistent allele coding across input GWAS. Required input columns (case-sensitive) are `SNP, A1, A2, BETA/Z/OR, SE, P, N, MAF`; some are conditional on `se.logit` / `OLS`. Silent failures are almost always column-name mismatches (e.g. `EA`/`NEA` instead of `A1`/`A2` -> 0% SNPs retained) or a missing `N` column -> SNPs dropped. Always run `head(read.table(file, header=TRUE, nrow=2))` per input before the `sumstats()` call.

```r
# Prepare per-SNP betas and SEs across all input GWAS
ss <- sumstats(
    files = c('raw/trait1.txt', 'raw/trait2.txt', 'raw/trait3.txt'),
    ref = 'reference.1000G.maf.0.005.txt',
    trait.names = trait_names,
    se.logit = c(TRUE, TRUE, FALSE),     # TRUE if trait is logistic-scale (case-control); FALSE if continuous
    OLS = c(FALSE, FALSE, TRUE),
    linprob = c(FALSE, FALSE, FALSE),
    N = N,
    info.filter = 0.9,                   # standard imputation INFO threshold
    maf.filter = 0.01
)

saveRDS(ldsc_results, 'ldsc_results.rds'); saveRDS(ss, 'sumstats_snps.rds')
```

```bash
# commonfactorGWAS DWLS + ML cross-check + Q_SNP classification; last argument = cores
Rscript scripts/commonfactor_gwas_qsnp.R ldsc_results.rds sumstats_snps.rds cfgwas_qsnp.tsv 8
```

Output: the leading `sumstats()` columns (`SNP`, `CHR`, `BP`, `MAF`, `A1`, `A2` ...), then `i`, `lhs`, `op`, `rhs`, `est`, `se_c`, `Z_Estimate`, `Pval_Estimate` (factor effect), `Q` / `Q_df` / `Q_pval` (`Q_pval` IS the per-SNP Q_SNP test), `fail`, `warning`, and the script's `factor_sig`, `qsnp_sig`, `factor_only_dwls`, `Q_pval_ML`, `factor_only`. Cluster and `MPI=` settings: `references/runtime-and-python-tools.md`.

The "factor-only" subset (factor-significant AND Q_SNP non-significant) is the publication-grade set of common-factor SNPs.

**Cross-check Q_pval with `estimation = 'ML'`; DWLS Q_pval is provisional.** On two synthetic panels (flat SE; MAF/N-driven SE) built without `ldsc()`, `commonfactorGWAS(estimation = 'DWLS')` returned Q_pval 0.91-0.96 for both 5 planted heterogeneous SNPs and 5 planted factor SNPs (no separation), while `'ML'` on the identical input separated them (heterogeneous 6e-83 to 6e-11, factor 0.28-0.48; GenomicSEM 0.0.5 + lavaan 0.6.19, re-audit Inputs 2 and 8). Whether that is a DWLS weakness or an artifact of synthetic V without per-SNP N is unsettled. So: when V/N come from anything other than a genuine `ldsc()` + `sumstats()` run, treat DWLS Q_pval as provisional; and even on real data, do not call a SNP factor-only on a DWLS Q_pval alone.

The script runs that ML pass itself and sets `factor_only` only when DWLS and ML Q_pval both clear the threshold.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `usermodel()`/`commonfactorGWAS()`/`userGWAS()` error `object 'ReorderModel'` or `'ReorderModelnoSNP' not found` (both DWLS and ML) | lavaan >=0.7.0 requires `ordered=FALSE` for DWLS on continuous data; GenomicSEM 0.0.5 never supplies it (see Version Compatibility) | Pin lavaan to 0.6.19 (see Tool Installation); confirmed working under both estimators |
| `commonfactor()` reports "failed to converge" but `traceback()` shows an underlying `object '...Results' not found` error | Same lavaan-version incompatibility above, mislabeled by GenomicSEM's tryCatch as non-convergence | Not a model-specification problem -- do not re-specify the model; pin lavaan to 0.6.19 |
| `commonfactor()` complains "S not positive definite" | Genetic correlations near +/-1 among inputs | Drop redundant traits; verify rg < 0.95 pairwise |
| `ldsc()` fails with "category not found" | Wrong LD score column names (legacy format) | Use Python 3 LDSC fork; download `eur_w_ld_chr/` from alkesgroup |
| `lavaan` says "model not identified" | Too few traits for too many parameters | Need >= 3 traits per factor; constrain factor variance to 1 |
| MTAG `MaxFDR` not in log | Older MTAG version (< 1.0.7) | Update MTAG; MaxFDR reporting added late 2019 |
| `usermodel()` slow or fails | Complex syntax + many traits | Simplify model; estimate with `estimation = 'DWLS'`, not `'ML'`, when V is informative |
| Sumstats output has zero overlap with reference | Allele coding mismatch in `sumstats()` | Check `se.logit` and `OLS` settings per trait; align A1/A2 |
| GenomicSEM and TwoSampleMR give different rg | TwoSampleMR uses bivariate LDSC; GenomicSEM uses the same S | Match the underlying LDSC reference panel and weights |

## Tool Installation

```r
# GenomicSEM is GitHub-only
remotes::install_github('GenomicSEM/GenomicSEM')

# Pin lavaan to 0.6.19 -- lavaan >=0.7.0 breaks usermodel()/commonfactorGWAS()/userGWAS()
# (see Version Compatibility). install_github above may pull a newer lavaan as a dependency;
# always run this line after it, and re-run it if packageVersion('lavaan') drifts to >=0.7.
remotes::install_version('lavaan', version = '0.6-19')

# Dependencies
install.packages(c('Matrix', 'gdata'))

# Optional companions
remotes::install_github('MRCIEU/TwoSampleMR')  # for downstream MR using factor GWAS as exposure
```

LDSC and MTAG Python installs: `references/runtime-and-python-tools.md`. Pre-downloaded reference files: `eur_w_ld_chr/`, `baselineLD_v2.2.*`, `w_hm3.snplist`, and 1000G allele-frequency files are hosted at `alkesgroup.broadinstitute.org/LDSCORE/`.

## References

- Grotzinger AD et al 2019 Nat Hum Behav 3:513 (GenomicSEM, common-factor GWAS, Q_SNP)
- Grotzinger AD et al 2022 Nat Genet 54:548 (Stratified GenomicSEM)
- Turley P et al 2018 Nat Genet 50:229 (MTAG; MaxFDR)
- Bulik-Sullivan B et al 2015 Nat Genet 47:291 (LDSC for genetic covariance)
- Bulik-Sullivan B et al 2015 Nat Genet 47:1236 (bivariate LDSC, sample overlap)
- Rosseel Y 2012 J Stat Softw 48:1-36 (lavaan package)
- Hu LT & Bentler PM 1999 Struct Equ Model 6:1 (CFI / RMSEA cutoffs)
- Asparouhov T & Muthen B 2009 Struct Equ Model 16:397 (ESEM framework)
- Finucane HK et al 2015 Nat Genet 47:1228 (S-LDSC, foundation for stratified GenomicSEM)
- Gazal S et al 2017 Nat Genet 49:1421 (baseline-LD annotations)
- Demange PA et al 2021 Nat Genet 53:35 (GenomicSEM GWAS-by-subtraction for noncognitive skills; Q_SNP in practice)
- Karlsson Linner R, Mallard TT et al 2021 Nat Neurosci 24:1367 (multivariate externalizing GWAS via GenomicSEM)
- de la Fuente J et al 2021 Nat Hum Behav 5:49 (GenomicSEM for cognitive g factor)
- Skrivankova VW et al 2021 JAMA 326:1614 (STROBE-MR; relevant when downstream MR uses factor GWAS)

## Related Skills

- causal-genomics/mendelian-randomization - Use the Q_SNP-clean "factor-only" subset of factor-GWAS effect sizes (see Common-Factor GWAS with Q_SNP) as MR exposure
- causal-genomics/genetic-correlation - Bivariate LDSC produces the off-diagonals of the S matrix; GenomicSEM is the multi-trait extension
- causal-genomics/heritability-partitioning - LDSC and S-LDSC foundations for stratified GenomicSEM
- causal-genomics/colocalization-analysis - Cross-trait colocalization at common-factor loci
- causal-genomics/pleiotropy-detection - Q_SNP is a per-SNP pleiotropy diagnostic; sibling concept
- causal-genomics/fine-mapping - Resolve factor-significant loci to credible sets
- causal-genomics/mediation-analysis - SEM mediation paths overlap with `usermodel` indirect effects
- causal-genomics/transcriptome-wide-association - TWAS on factor sumstats from `commonfactorGWAS`
- population-genetics/association-testing - GWAS sumstats are the input format
