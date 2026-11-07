# Per-Method Failure Modes
### Cross-trait LDSC intercept misread as bias

**Trigger:** Reader (collaborator, reviewer) sees a non-zero LDSC intercept and reports the rg as "biased by sample overlap".

**Mechanism:** The bivariate LDSC model partitions covariance between traits into a slope (rg-driven, scales with LD score) and an intercept (overlap-driven, constant in LD score). The slope is what carries rg, and it remains unbiased regardless of intercept value (Bulik-Sullivan 2015 Nat Genet 47:1236, Methods). Confusing the intercept with bias on the slope is a routine misinterpretation.

**Symptom:** Reviewer comment requesting "correction for sample overlap" when LDSC was already used; collaborator suggesting switching to HDL because intercept is non-zero.

**Fix:** Report rg with SE and the intercept as a separate statistic; cite Bulik-Sullivan 2015 Methods explicitly; do NOT switch to HDL (which is the wrong direction since HDL is the method that breaks under overlap, not LDSC).

### HDL bias with sample overlap

**Trigger:** Running HDL on two GWAS that share > 5% of individuals (e.g. two UKB-derived GWAS, two MVP-derived GWAS, GWAS reusing controls).

**Mechanism:** HDL likelihood assumes independent trait residuals after marginalizing genetics. Overlap induces non-zero residual correlation; the likelihood is misspecified and the estimate is pulled toward phenotypic correlation in the shared subset (Ning 2020 Nat Genet 52:859 Supplement).

**Symptom:** HDL rg differs substantially from cross-trait LDSC rg; HDL CI is narrower than expected from N alone; LDSC intercept (which absorbs overlap) is materially non-zero.

**Fix:** Compute LDSC intercept first as an overlap diagnostic; if non-zero, switch to LDSC as primary. HDL remains valid only when the two GWAS draw from non-overlapping cohorts (verify by cohort identifier, not just by file source).

### Non-EUR ancestry mismatch with EUR LD scores

**Trigger:** Running LDSC on a non-EUR GWAS (or admixed sample) with the default `eur_w_ld_chr/` reference.

**Mechanism:** LD scores are ancestry-specific; mean LD per SNP differs across populations and bivariate moments use the wrong null. h2 and rg estimates are systematically biased; intercept can be inflated.

**Symptom:** LDSC ratio is unusually high (>0.2); per-chromosome estimates wildly heterogeneous; total h2 mismatches independent estimates from the same cohort.

**Fix:** Use ancestry-matched LD scores from alkesgroup (eas_w_ld_chr, afr_w_ld_chr) or compute custom LD scores from in-sample LD reference. For trans-ancestry rg, switch to Popcorn.

### Global rg masks local rg variation

**Trigger:** Two traits with biologically plausible shared etiology return global rg ~ 0 in cross-trait LDSC.

**Mechanism:** Global rg averages over the genome; loci with positive local rg can cancel against loci with negative local rg, particularly for traits with antagonistic pleiotropy (e.g. autoimmune vs infectious-disease susceptibility), or when shared etiology is confined to a small fraction of the genome.

**Symptom:** Well-powered (mean chi-square >> 1.02) global LDSC rg near zero with wide CI overlapping zero, while domain biology, prior co-occurrence studies, or shared-pathway analyses strongly suggest shared etiology.

**Fix:** Run LAVA over the standard ~2495 LDetect loci; per-locus Bonferroni-significant rg at any locus is evidence of localized shared etiology. Annotate hit loci with overlapping GWAS catalog signals and pathway/tissue enrichment.

### Cross-population rg below 1 even at causal level

**Trigger:** Computing rg between same-trait GWAS in two ancestries (e.g. EUR T2D vs EAS T2D).

**Mechanism:** Even when the trait has the same biological definition, causal variant identity and effect sizes can differ across populations due to gene-environment interaction, allele-frequency divergence, and population-specific epistasis. Popcorn 2016 demonstrated that rg(cross-pop) < 1 is common and biologically real, not a methodological artifact.

**Symptom:** Trans-ancestry rg point estimate around 0.6-0.9 with CI excluding 1 for a trait expected to be "the same disease".

**Fix:** Use Popcorn (within-pop LDSC is invalid for cross-pop rg); interpret rg(cross-pop) < 1 as quantifying population-specific architecture rather than as bias; report rho_ge (effect correlation) and rho_gi (impact correlation) separately.

### Same-Trait Cross-Cohort rg as Consistency Check

**Use case:** Meta-analysis QC -- two same-trait GWAS (e.g. UKB IBD vs FinnGen IBD) yield rg < 1 with CI excluding 1. This is distinct from the cross-population analog above (handled by Popcorn); here both cohorts are same-ancestry but different studies.

**Interpretation:** (a) population-substructure differences, (b) phenotype-definition heterogeneity (e.g. different ICD coding, self-report vs registry), (c) genuine biology (founder effects in isolates like FinnGen).

**Decision rule:** rg ~ 0.9-1.0 with CI overlapping 1 -> consistent enough to meta-analyze; rg ~ 0.7-0.9 -> moderate heterogeneity, consider sensitivity meta with random effects; rg < 0.7 -> re-examine phenotype definitions before meta-analyzing.

### Low chi-square mean (underpowered GWAS)

**Trigger:** Mean chi-square in either input GWAS is below 1.02 (heuristic LDSC threshold).

**Mechanism:** LDSC, HDL, and LAVA all depend on bivariate moments of Z-scores against LD score; weak signal means the slope is dominated by noise.

**Symptom:** LDSC rg SE > 0.2; intercept estimates fluctuate across chromosome; HDL convergence warnings; LAVA returns p > 0.05 at most loci.

**Fix:** Meta-analyze contributing cohorts to push effective N to ~50k or higher before running rg; if meta-analysis is not feasible, report rg as exploratory; do NOT switch methods to "salvage" power, the problem is upstream of method choice.

### LAVA univariate filter ignored

**Trigger:** Reporting LAVA bivariate rg at a locus where the univariate local h2 is non-significant in one or both traits.

**Mechanism:** LAVA's bivariate test is only valid at loci with detectable local heritability in BOTH traits. Without local h2 signal in at least one trait, the bivariate test is unidentified and can return spurious significant rg.

**Symptom:** LAVA bivariate p < 0.05 at loci where univariate h2 p > 0.05 for one trait; rg estimates near +/- 1 (boundary cases).

**Fix:** Filter loci on `univ.p < 0.05 / N_loci` (Bonferroni for ~2495 loci) in BOTH traits before running `run.bivar()`; report only at filtered loci. This is the documented LAVA workflow in Werme 2022 Supplement and the GitHub vignette.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| |rg| > 0.7 | Operational high correlation | Suggests strong shared genetic architecture; near-universal in pairs like MDD-anxiety or LDL-CHD |
| |rg| 0.3-0.7 | Operational moderate correlation | Common in psychiatric / cardiometabolic trait pairs; routine to flag for joint analysis |
| |rg| < 0.3 | Operational low correlation | Globally weak; may still harbor biologically meaningful local rg via LAVA |
| rg SE < 0.05 | Operational reliable estimate | Above this SE, point estimate is uncertain to 1 decimal place |
| HDL sample overlap < 5% | Ning 2020 Nat Genet (Supplement) | Above this, HDL likelihood is misspecified and biased |
| LDSC mean chi-square > 1.02 | LDSC documentation (Bulik-Sullivan 2015 tutorial) | Below this, LD-score regression is severely underpowered for h2 / rg |
| LAVA local p < 0.05 / N_loci | Werme 2022 Nat Genet 54:274 | Bonferroni for ~2495 LDetect loci; standard genome-wide local-rg correction |
| Popcorn rho_ge CI excludes 1 | Brown 2016 AJHG 99:76 | Evidence of population-specific causal architecture |
| LCV gcp != 0 (two-sided p < 0.05) | O'Connor & Price 2018 Nat Genet 50:1728 | Directional evidence of (partial) causation given non-zero rg |
| LDSC ratio < 0.2 | Bulik-Sullivan 2015 Methods | High ratio (intercept / chi-square - 1) indicates population stratification or model misfit |
| MR + rg sensitivity trigger | Operational | |rg| > 0.3 with a significant IVW estimate REQUIRES CHP-aware sensitivity (CAUSE / LHC-MR) |
| Conditional-rg LAVA covariate set | Werme 2022 | Up to 4 conditioning traits per `run.pcor()` call before identification fails |
