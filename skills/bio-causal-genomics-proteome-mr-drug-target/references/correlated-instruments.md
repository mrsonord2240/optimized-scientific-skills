# Correlated cis-pQTL Instruments

Moved verbatim from SKILL.md (one pointer edited: `dat` comes from `examples/cis_pqtl_mr.R`). Read when several cis-pQTLs in the window are correlated (r2 0.1 to 0.7), for the `mr_ivw` namespace and sign-alignment caveats, or for the robust and Patel 2023 alternatives.

## Cis-IVW with Correlated Instruments

**Goal:** When several cis-pQTLs in the window are correlated (r2 between 0.1 and 0.7), use the generalized IVW that takes the LD matrix as an explicit parameter.

**Approach:** Compute or load the LD matrix in the cis-window from an ancestry-matched plink reference; build the input with `MendelianRandomization::mr_input(..., correlation = ld_matrix)` then call `MendelianRandomization::mr_ivw(mr_obj, model = 'default')`.

```r
library(MendelianRandomization); library(ieugwasr)

# Harmonised SNPs only, so the matrix rows match the beta vectors (dat comes from `examples/cis_pqtl_mr.R`)
ld <- ld_matrix(dat$SNP, bfile = '1kg_EUR/EUR', plink_bin = genetics.binaRies::get_plink_binary())

mr_obj <- mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                   by = dat$beta.outcome, byse = dat$se.outcome,
                   correlation = ld)

result_correl <- MendelianRandomization::mr_ivw(mr_obj, model = 'default')
```

Numerical caveat: when any pair of cis-pQTLs has r2 ~ 1 (e.g. perfect proxies), the LD matrix is rank-deficient and the SE explodes. Pre-prune at r2 < 0.95.

**API caveat for `MendelianRandomization::mr_ivw`:**

- **Always namespace it.** `TwoSampleMR::mr_ivw` is a plain function with a different signature. When TwoSampleMR is attached after MendelianRandomization (the normal order in this Skill, together with coloc), the bare `mr_ivw(mr_obj, model = 'default')` resolves to the TwoSampleMR one and fails with `unused arguments`. Checked on MendelianRandomization 0.10.0 / TwoSampleMR 0.7.9: `environmentName(environment(mr_ivw))` returned `TwoSampleMR` after `library(MendelianRandomization); library(TwoSampleMR)`.
- Correlation between cis-pQTLs is supplied at MRInput construction via `correlation=` (`corr=` only works through R partial-matching; use the full name). `mr_ivw()` reads that slot itself. Checked on 0.10.0: the estimate and SE are identical with `correl = TRUE`, `correl = FALSE` and no `correl` when the matrix is present, and `correl = TRUE` with no matrix errors. Do not pass `correl`.
- **Sign alignment.** `ld_matrix()` returns signed r relative to the 1000G major alleles (ieugwasr docs), and `mr_ivw()` prints a warning that correlations must be relative to the same effect alleles as the estimates. `with_alleles = TRUE` (the default) appends the reference alleles to the row names so the signs can be checked against `dat$effect_allele.exposure`; flip rows and columns whose reference allele differs, then strip the allele suffix before `mr_input()`.

### Robust / Penalized cis-IVW as a Correlated-Instrument Sensitivity Estimator

**Goal:** Provide a robust sensitivity estimate when cis-pQTLs are correlated and a minority may be outliers.

**Approach:** `MendelianRandomization::mr_ivw(robust=TRUE, penalized=TRUE)` applies Burgess's robust-regression + penalized-weights IVW, down-weighting heterogeneous/outlying instruments. A distinct, more modern option is Patel, Gill, Newcombe, Burgess 2023 *Biometrics* 79:3458-3471, which reduces the dimension of correlated cis-variants in a single gene region via factor analysis and applies weak-factor-robust conditional inference; it exploits the within-region genetic-correlation (LD) structure rather than avoiding it, and is implemented separately from `mr_ivw`.

```r
library(MendelianRandomization)

mr_obj <- mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                   by = dat$beta.outcome, byse = dat$se.outcome,
                   correlation = ld)
robust_res <- MendelianRandomization::mr_ivw(mr_obj, model = 'default', robust = TRUE, penalized = TRUE)
```

Decision rule: standard cis-IVW when post-clumped LD r2 < 0.1; correlated-IV cis-IVW (Burgess 2017 Genet Epidemiol) when r2 between 0.1 and 0.7 AND ancestry-matched LD matrix is trustworthy; robust/penalized IVW (above) as a sensitivity check when a minority of instruments may be outliers; Patel 2023 conditional cis-MR when instruments are highly correlated and weak (factor-analysis + conditional inference), and preferred when the LD reference is ancestry-mismatched. Benchmarks for these correlated-IV estimators are still evolving; report more than one as sensitivity when feasible.

