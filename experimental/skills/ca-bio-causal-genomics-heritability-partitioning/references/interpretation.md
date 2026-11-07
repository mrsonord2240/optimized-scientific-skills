# LDSC Intercept Interpretation (Postdoc Nuance)
The LDSC intercept is widely misinterpreted as a "confounding score". The correct interpretation:

- Intercept = 1 indicates no inflation from population structure, cryptic relatedness, or sample overlap (idealised)
- Intercept > 1 indicates SOME source of inflation, BUT polygenic background can elevate the intercept too: at very high N, polygenic signal can lift the intercept modestly without stratification
- The **ratio** statistic `ratio = (intercept - 1) / (mean_chi2 - 1)` is the fraction of inflation attributable to non-polygenic sources; a ratio of 0 means all inflation is polygenic, a ratio near 1 means most of it is stratification or overlap
- Bulik-Sullivan 2015 recommends interpreting intercept jointly with mean chi-square; intercept of 1.05 is innocuous if mean chi-square is 1.5 (ratio = 0.1) but worrying if mean chi-square is 1.05 (ratio = 1.0)

**Operational rule:** Always report intercept, mean chi-square, and ratio together. Do not interpret intercept in isolation. For sample-overlap diagnosis between two GWAS, use bivariate LDSC intercept, not univariate.

| Intercept | Ratio | Interpretation |
|-----------|-------|----------------|
| ~1.0 | ~0 | No inflation; h2 trustworthy |
| 1.0 - 1.1 | < 0.2 | Mostly polygenic; h2 trustworthy |
| 1.1 - 1.3 | 0.2 - 0.5 | Mild inflation; investigate population structure / sample overlap |
| 1.3 - 1.5 | 0.5 - 0.8 | Substantial inflation; report ratio jointly; consider re-genotype-QC |
| > 1.5 | -- | Re-run after PC adjustment or genomic control; do not interpret h2 |

**Intercept > 1.5 troubleshooting ladder** (work in order; stop when source is found): (a) per-cohort PC adjustment was insufficient; refit GWAS with more PCs (10-20) or per-cohort separately, (b) cryptic relatedness in the GWAS cohort -- run `king --related` and remove pairs with kinship > 0.05 (or 0.0884 for second-degree), (c) case-control matching imbalance -- check case/control PCs separately, (d) sample-overlap with one of the contributing cohorts (especially in meta-analysis) -- check bivariate intercepts pairwise, (e) if biobank-internal, recompute the GWAS with sample-level relatedness exclusion before LDSC.

## LDSC vs LDAK Reconciliation

LDSC (GCTA model) and LDAK SumHer (LDAK-Thin model) make different assumptions about how per-SNP heritability scales with LD and MAF:

- **GCTA model (LDSC default):** per-SNP h2 inversely proportional to local LD score; high-LD SNPs tag many causal variants
- **LDAK-Thin (Speed 2017 Nat Genet 49:986-992 introduced the LDAK model; the Thin model was formalized in Speed 2020 Nat Genet 52:458):** per-SNP h2 weighted by MAF and inversely by LD with empirical exponents; less weight on common high-LD SNPs

These give systematically different functional enrichment estimates. Speed 2019 reported that LDAK-Thin often gives **lower** conserved-region enrichment than baseline LDSC; conversely LDSC can over-attribute h2 to coding/conserved regions because the GCTA prior couples LD to causality. Gazal 2019 Nat Genet 51:1202-1204 argues the baseline-LD S-LDSC model is better-calibrated (higher model likelihood) and that LDAK/SumHer's lower functional-enrichment estimates are downward-biased; the LDAK developers dispute this (Speed 2020 Nat Genet 52:458), so report both models and flag the model-dependence.

**Operational rule:** Whenever functional enrichment is the primary claim (e.g. "h2 is enriched in tissue T by N-fold"), report enrichment from BOTH LDSC and LDAK. Flag the model assumption. If LDSC and LDAK disagree by > 2x (e.g. LDSC conserved-region enrichment 25x vs LDAK 10x; both are model-internally consistent and the data alone cannot pick one), treat the claim as model-dependent and cite Gazal 2019. Report LDSC primary + LDAK confirmatory, emphasize directional agreement over magnitude, and never pick the model that gives the desired answer. For non-enrichment claims (total h2, rg between two traits), the model dependence is smaller and LDSC alone is acceptable.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| mean chi-square > 1.02 | LDSC wiki / Bulik-Sullivan 2015 | Below this, LDSC h2 estimate has huge SE; need N proportional to 1 / h2 |
| Total h2 must lie in [0, 1] on both observed and liability scale | Definition of h2_SNP | A point estimate outside [0, 1] is never a face-value result. If the 95% CI (estimate +/- 1.96 SE) still excludes [0, 1] (e.g. liability h2 1.43 (0.10)), it is an artifact: check `--samp-prev`/`--pop-prev` assignment, munged N and the LD reference, then report it as a data or model failure. If the CI overlaps the bound, report it as noise around the boundary |
| h2 SE < 0.02 | LDSC convention | Below this SE, h2 estimate is interpretable; above, treat as exploratory |
| h2 SE >= 0.02 OR mean chi-square < 1.02 | LDSC convention | Estimate unreliable; N >= 50k is typical noise floor for h2 ~ 0.1 (scales as ~1/h2) |
| LDSC intercept / ratio bins | Bulik-Sullivan 2015 | See the intercept/ratio table under "LDSC Intercept Interpretation" |
| Stratified LDSC enrichment p < 0.05 / N_annot | Finucane 2015 | Bonferroni across annotations in the baseline-LD model |
| S-LDSC cell-type p < 2.5e-4 | Finucane 2018; ~200 tissues | Bonferroni for tissue prioritization |
| HESS h2 per locus needs >= 1000 SNPs | Shi 2017 AJHG 101:737 | Quadratic form unstable below this density |
| HDL sample overlap < 5% | Ning 2020 Nat Genet 52:859 | Above 5%, HDL likelihood is biased |
| LDAK tagging file build match | Speed 2019 | hg19 vs hg38 tagging files non-interchangeable |
| Annotation > 0.5% of genome | Finucane 2015 | Smaller categories underpowered for tau estimation |
| Effective N > 5000 per population (Popcorn) | Brown 2016 AJHG 99:76 | Below this, trans-ancestry rg has very wide CI |
