# LDSC and S-LDSC failure modes (reference for heritability-partitioning SKILL.md)

### LDSC intercept misinterpretation

**Trigger:** Reporting intercept ~1.05 as "evidence of confounding".

**Mechanism:** Intercept absorbs mean chi-square inflation from any non-polygenic source PLUS some polygenic contribution at very high N. In isolation, intercept above 1 does not imply confounding.

**Symptom:** Methods section claims population stratification based solely on intercept value; reviewers flag the omission of ratio statistic.

**Fix:** Always report intercept, mean chi-square, ratio = (intercept - 1) / (mean_chi2 - 1), and h2 jointly. Interpret ratio < 0.2 as "mostly polygenic, h2 trustworthy"; ratio > 0.3 as "investigate population structure / overlap before claiming h2".

### LDSC with non-EUR ancestry and EUR LD scores

**Trigger:** Applying default EUR LD scores from `alkesgroup.broadinstitute.org/LDSCORE/eur_w_ld_chr/` to an EAS, AFR, or AMR GWAS.

**Mechanism:** LD-score regression assumes the LD-score covariate matches the GWAS population's LD structure. Cross-ancestry application produces biased h2 (typically underestimates) and inflated intercept.

**Symptom:** h2 estimate < 0.05 despite trait being known-heritable from twin / family studies; intercept > 1.2 with non-polygenic mean chi-square; ratio > 0.5.

**Fix:** Use ancestry-matched LD scores (EAS, AFR, AMR available at the same Alkes group URL). If multi-ancestry meta-analysis, use Popcorn or trans-ancestry MAMA framework rather than LDSC on the combined sumstats.

### Stratified LDSC with collinear annotations

**Trigger:** Adding a custom annotation that overlaps heavily with an existing baseline category (e.g. "active promoter" against "promoter").

**Mechanism:** Per-annotation tau coefficients are estimated jointly via multivariable regression; collinearity inflates per-tau SE and can flip the sign of marginal effect.

**Symptom:** Custom annotation tau has very large SE, p-value > 0.5; baseline categories that were significant become non-significant.

**Fix:** Test annotations marginal to the baseline by including baseline-LD_v2.2 plus the new annotation only; never test multiple highly correlated annotations jointly; report VIF of annotation matrix; use the joint enrichment of {baseline + new} category not per-tau.
