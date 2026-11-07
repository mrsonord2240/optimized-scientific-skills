# Model Fit Diagnostics
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
