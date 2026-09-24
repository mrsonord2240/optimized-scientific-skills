# Advanced model reference (reference for genomic-sem SKILL.md)

ESEM, `userGWAS()` custom path models, and higher-order / bifactor / p-factor models. Read it when `SKILL.md` "Reference Files" points here. The section names below refer to `SKILL.md`.

## ESEM (Exploratory Factor Structure)

When the factor structure is unknown, fit `usermodel()` with all loadings free across all factors, then apply a rotation post-fit. Rotation choices: **geomin oblique** (default; allows factor correlation), **target rotation** (Browne 2001 Multivariate Behav Res 36:111; uses a hypothesized loading template), **quartimin** (orthogonal; assumes factors are uncorrelated).

```r
model_esem <- '
    F1 =~ NA*trait1 + trait2 + trait3 + trait4
    F2 =~ NA*trait1 + trait2 + trait3 + trait4
    F1 ~~ 1*F1
    F2 ~~ 1*F2
    F1 ~~ F2
'
esem_fit <- usermodel(covstruc = ldsc_results, model = model_esem, estimation = 'DWLS')
# Rotate post-fit via GPArotation::GPForth/GPFoblq or lavaan::rotate()
```

**Decision:** ESEM for K-factor exploration when structure is unknown; CFA via `usermodel()` once a structure is confirmed. Report rotation sensitivity (geomin vs target vs quartimin) and treat as exploratory.

**Cross-loadings.** Cross-loadings (one trait loads on > 1 factor) are common in psychiatric and behavioral GWAS. Brown 2015 *Confirmatory Factor Analysis for Applied Research* recommends allowing cross-loadings first and using modification indices to guide simplification. Allow a cross-loading when constraining residual variance otherwise forces a Heywood case. Constrain when CFI < 0.9 and modification indices instead suggest a correlated residual between two indicators (which is the more parsimonious fix).

## userGWAS for Custom Path Models

`userGWAS()` fits arbitrary lavaan-syntax SNP regressions and is the right tool when the SNP needs to be tested on multiple paths simultaneously (e.g. factor-mediated effect AND a direct effect on one indicator).

```r
# Test SNP -> F path + SNP -> trait1 direct path simultaneously
model <- '
    F =~ NA*trait1 + trait2 + trait3
    F ~~ 1*F
    F ~ SNP
    trait1 ~ SNP    # direct effect on trait1, partialed out of F
'
user_results <- userGWAS(covstruc = ldsc_results,
                         SNPs = ss,
                         estimation = 'DWLS',
                         model = model,
                         sub = c('F~SNP', 'trait1~SNP'),
                         parallel = TRUE,
                         cores = 8)
# Output (one data frame per SNP set; checked on 0.0.5): SNP, CHR, BP, MAF, A1, A2, lhs, op, rhs,
# free, label, est, SE, Z_Estimate, Pval_Estimate, chisq, chisq_df, chisq_pval, AIC, error, warning.
# There is NO Q_pval column: chisq / chisq_pval is the fit of the whole SNP-augmented model, so a
# small chisq_pval means unmodelled SNP paths remain (a third path may be needed).
```

## Higher-Order / Bifactor / p-Factor Models

Use case: psychiatric genetics p-factor (Caspi 2014 Clin Psychol Sci 2:119; Grotzinger 2022 Nat Genet 54:548 cross-disorder), cognitive g-factor (de la Fuente 2021 Nat Hum Behav 5:49).

Hierarchical template -- first-order factors load on a single second-order p-factor:

```r
model_pfactor <- '
    # First-order factors
    INT =~ NA*trait_anx + trait_dep + trait_neuro       # internalizing
    EXT =~ NA*trait_adhd + trait_alc + trait_subst       # externalizing
    THT =~ NA*trait_scz + trait_bp                       # thought-disorder
    # Second-order p-factor
    p =~ NA*INT + EXT + THT
    INT ~~ 1*INT
    EXT ~~ 1*EXT
    THT ~~ 1*THT
    p ~~ 1*p
'
```

The second-order p-factor follows the same rule as any factor: it needs >= 3 first-order factors to identify. With only 2 first-order factors, `usermodel()` returns chi-square ~ 0 and warns that the information matrix could not be inverted, so SEs are NaN (checked on 0.0.5 + lavaan 0.6.19, DWLS); drop the second-order factor or add a third first-order factor.

Bifactor alternative: `p =~` all traits directly, with `INT`/`EXT`/`THT` as orthogonal residual factors. Bifactor typically gives tighter CFI/RMSEA but the substantive interpretation of the residual factors is harder; bifactor is also prone to over-fitting at modest trait counts (Bonifay W & Cai L 2017 Multivariate Behav Res 52:465). Cite Grotzinger 2022 Nat Genet 54:548 and Karlsson Linner R, Mallard TT et al 2021 Nat Neurosci 24:1367 for the canonical psychiatric implementations.
