# Per-Method Failure Modes

## Per-Method Failure Modes

### MR-Egger NOME violation

**Trigger:** I^2_GX = (Q_GX - df) / Q_GX is below 0.9, indicating measurement-error attenuation of the Egger slope (NOME = "no measurement error" in the exposure GWAS effect sizes).

**Mechanism:** MR-Egger regresses outcome effects on exposure effects with a free intercept. Imprecise exposure effects (high beta.exposure SE relative to beta.exposure variability across instruments) introduce regression dilution that pulls the Egger slope toward the null and inflates the intercept.

**Symptom:** Egger slope much closer to zero than IVW, weighted median, and weighted mode estimates; large Egger SE.

**Fix:** Apply SIMEX correction (Bowden 2016 IJE 45:1961; Cook & Stefanski 1994 JASA 89:1314 SIMEX framework) using the `simex` package on the Egger regression, treating beta.exposure SE as measurement error. See examples/simex_egger_correction.R. Alternative: use MR-RAPS, which models the exposure-effect error explicitly via profile likelihood and does not suffer the NOME failure.

**Caveat:** SIMEX removes a known bias direction but its extrapolation step is noisy in small samples, so the corrected slope is not guaranteed to be closer to the truth than the naive one (in a 20-SNP audit run with I^2_GX 0.70 and planted effect 0.3, naive Egger 0.666 became SIMEX 0.854, while IVW 0.314 and contamination mixture 0.343 sat near truth). Report the naive and SIMEX-corrected slopes side by side with CIs and read them against IVW / median / mode; do not present the corrected value as strictly superior.

### MR-PRESSO majority-outlier breakdown

**Trigger:** More than 50% of instruments are pleiotropic (UHP), e.g. when instrument set was loosely selected (genome-wide significant but unfiltered).

**Mechanism:** MR-PRESSO's global RSS-out statistic and outlier detection both assume a majority-valid set; outliers are defined relative to that majority. With a pleiotropic majority, PRESSO removes the valid minority.

**Symptom:** PRESSO-corrected estimate is similar in magnitude (and sign) to the uncorrected estimate even after dropping nominally "outlier" SNPs; distortion-test p-value paradoxically non-significant; few or no outliers detected despite obvious global-test significance.

**Fix:** Do not trust PRESSO corrected estimate. Re-examine instrument selection (drop loose p-thresholds, prune LD harder); switch to CAUSE or LHC-MR; consider weighted-mode estimator which is plurality-valid rather than majority-valid.

### MR-PRESSO false negative under CHP

**Trigger:** Strong shared heritable confounder (high rg) producing CHP. Confirmed by significant LDSC rg or LCV gcp.

**Mechanism:** Correlated pleiotropy is a population-level mean shift in alpha conditional on gamma; it is not an outlier pattern. PRESSO's RSS-out distance is invariant under such a mean shift, so the global test is not powered against CHP.

**Symptom:** PRESSO global p > 0.05 (no detected pleiotropy) while a CHP-aware method (CAUSE, LHC-MR) returns a substantially different (often null) causal estimate.

**Fix:** When CHP is plausible, ALWAYS run CAUSE or LHC-MR in addition to PRESSO; do not rely on PRESSO global non-significance as evidence of no pleiotropy.

### MR-Egger underpowered with few SNPs

**Trigger:** Fewer than 10 instruments.

**Mechanism:** Egger's intercept variance is driven by the spread of beta.exposure across instruments; with few SNPs the intercept CI is so wide that even strongly pleiotropic data give non-significant intercepts.

**Symptom:** Non-significant Egger intercept p-value alongside obviously discordant IVW and weighted-median estimates.

**Fix:** Report intercept point estimate and CI rather than a binary "pleiotropy present / absent" verdict; do not use Egger as the only sensitivity method when SNP count is low; weight evidence toward weighted-median, weighted-mode, and CAUSE / LHC-MR.

### Steiger filter inverted by exposure measurement error (Hemani 2017)

**Trigger:** Exposure is imprecisely measured (lower heritability ascertained in the exposure GWAS) and outcome is well-measured.

**Mechanism:** Steiger compares r^2_GX vs r^2_GY per SNP. Measurement error in the exposure underestimates r^2_GX; well-measured outcome captures r^2_GY accurately. Per-SNP, the inequality can flip even when the true causal direction is exposure -> outcome.

**Symptom:** A large fraction of instruments fail Steiger (`steiger_dir == FALSE`) in a direction that conflicts with biological plausibility.

**Fix:** Interpret Steiger as one signal among many, not a hard gate; cross-check with bidirectional MR; verify exposure GWAS heritability and sample size; switch to LHC-MR which models both directions jointly and accounts for heritability.

### CAUSE underpowered with few significant SNPs

**Trigger:** Fewer than 100 genome-wide-significant instruments (p < 5e-8) after harmonization and LD pruning.

**Mechanism:** CAUSE fits a Bayesian mixture model over a shared-factor (CHP) component, a shared-causal component, and a null component. Posterior identification of the mixture weights requires substantial signal across many SNPs.

**Symptom:** CAUSE delta_ELPD CI crosses zero; Pareto-k diagnostic flags unstable points; posterior intervals on q (CHP fraction) span [0, 1].

**Fix:** Use LCV gcp for genome-wide directional inference (does not require many significant SNPs); use LHC-MR if heritability and sumstats are available; or report CAUSE alongside an explicit caveat about its underpowered regime.

### LCV gcp under non-Gaussian effect distributions

**Trigger:** Highly polygenic trait with substantial sparsity in true effects (mixture of large-effect and zero-effect loci).

**Mechanism:** LCV assumes a bivariate normal model for effect sizes after LDSC adjustment. Sparse architectures (e.g. immune traits with HLA dominance) violate this and bias gcp estimates.

**Symptom:** LCV gcp point estimate appears extreme but heritability LDSC z-scores are modest; partitioned heritability shows extreme HLA enrichment.

**Fix:** Exclude HLA region from LDSC inputs; complement with CAUSE / LHC-MR; report gcp with awareness of the polygenicity caveat.
