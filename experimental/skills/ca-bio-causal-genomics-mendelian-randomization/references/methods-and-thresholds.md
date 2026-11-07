# Per-Method Failure Modes
### IVW under directional pleiotropy

**Trigger:** Several SNPs affect the outcome through pathways not via the exposure, in a consistent direction.

**Mechanism:** IVW is a weighted regression through origin; non-zero mean pleiotropy shifts the slope.

**Symptom:** Egger intercept p < 0.05 with non-zero estimate; IVW differs from weighted median; MR-PRESSO global test p < 0.05.

**Fix:** Use Egger (if `I^2_GX >= 0.9` -- otherwise SIMEX-correct; code and caveats under "NOME violation invalidating Egger" below); cross-check with weighted median, MR-PRESSO, and CAUSE; report IVW only as one of a panel, never alone. The `MendelianRandomization::mr_egger()` function accepts `distribution='normal'` and reports the `I.sq` (I^2_GX) NOME diagnostic but applies no NOME/SIMEX correction to the estimate itself and does NOT expose a SIMEX wrapper.

### Weak-instrument bias direction

**Trigger:** Mean per-instrument F-statistic < 10, or several individual F < 10.

**Mechanism:** Weak IVs amplify finite-sample correlation between IV-X and IV-Y errors; bias direction depends on overlap regime (see table above).

**Symptom:** Estimates shift markedly when removing the weakest instruments; one-sample MR estimates much larger than two-sample.

**Fix:** Compute F per instrument from the EXPOSURE GWAS, not the outcome; exclude F < 10; use MR-RAPS (handles weak IVs by design); for two-sample, also report unweighted IVW (less weak-IV-bias-inflated than weighted in some regimes).

### Winner's curse at P~5e-8

**Trigger:** Discovery GWAS is the source of both instrument selection and effect-size estimates.

**Symptom:** MR effect shrinks substantially when using effect sizes from an independent replication GWAS.

**Mechanism and fix:** regression toward the mean inflates `beta_X`; three-sample design, sumstats-only corrections (MRlap, MR-SimSS, RIVW) and the Jiang 2023 inflation magnitudes: `references/mrlap-overlap-correction.md`.

### NOME violation invalidating Egger

**Trigger:** Running MR-Egger with `I^2_GX < 0.9`.

**Mechanism:** Egger assumes NO Measurement Error in exposure effect sizes (NOME); when violated, Egger slope is attenuated toward null with reciprocal bias on the intercept.

**Symptom:** `mr_pleiotropy_test()` Egger estimate disagrees with weighted median in magnitude but agrees in direction; `Isq()` function returns <0.9.

**Fix:** Compute `Isq(beta_X, se_X)` (Bowden 2016 IJE 45:1961); if <0.9, apply SIMEX correction via `simex` package or report Egger as exploratory only. The MendelianRandomization package's `mr_egger()` reports the `I.sq` (I^2_GX) NOME diagnostic but does not itself apply a NOME/SIMEX correction; SIMEX must be run separately.

Runnable SIMEX code (with the precomputed-weights workaround for the simex 1.8 crash): `references/simex-egger-nome.md`.

### Steiger filter false flag under unmeasured confounding

Steiger is a heuristic; it can flag "reverse causation" wrongly when an unmeasured confounder (e.g. SES) sits upstream of both traits (Lutz 2022). Cross-check direction with bidirectional MR, or use LCV / LHC-MR for confounder-rich traits. Trigger, mechanism and fix: `references/bidirectional-steiger.md`.

### Palindromic SNP harmonization

**Trigger:** SNPs with alleles A/T or C/G near MAF 0.5.

**Mechanism:** Strand orientation is ambiguous for palindromic SNPs when allele frequencies are intermediate; flipping introduces sign errors that look like pleiotropy.

**Symptom:** `harmonise_data()` reports many palindromic SNPs dropped; remaining SNPs show heterogeneity from a handful.

**Fix:** Default `action = 2` (infer from allele frequencies) drops MAF~0.5 palindromes; `action = 3` drops ALL palindromes (most conservative); never use `action = 1` (assumes forward strand) unless both GWAS are guaranteed to use the same strand convention. Document choice in methods.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| F-statistic > 10 per instrument | Staiger & Stock 1997 (linear IV) | Heuristic; debated (Burgess 2011 IJE; Zhao 2020 argues 10 is too low for one-sample) |
| Conditional F > 10 per exposure (MVMR) | Sanderson 2019 IJE 48:713 | Total F can be high while conditional F low; per-exposure F is what matters |
| `I^2_GX >= 0.9` for Egger | Bowden 2016 IJE 45:1961 | NOME assumption; below this, SIMEX correction required |
| CAUSE >= 100 significant SNPs | Morrison 2020 Nat Genet 52:740 | Mixture model needs signal density for shared-factor estimation |
| Egger >= 10 instruments | Bowden 2015 IJE 44:512 | Power for slope test in weighted regression |
| Sample-overlap z-score correlation | Burgess 2016 Genet Epidemiol 40:597 | Use LDSC bivariate intercept as proxy; correct IVW SE accordingly |
| Clumping r2 < 0.001, 10 Mb window | TwoSampleMR default; matches GWAS LD norms | Polygenic MR; cis-MR uses r2 < 0.1 within window |
| Steiger p < 0.05 | Hemani 2017 PLoS Genet 13:e1007081 | Heuristic; subject to confounder caveat (Lutz 2022) |
| MR-PRESSO NbDistribution | 1000 exploratory; >= 5000 publication; >= 10000 stringent | Verbanck 2018 Nat Genet 50:693; precision of global p-value scales with NbDistribution |
| STROBE-MR all 20 items | Skrivankova 2021 JAMA 326:1614; BMJ 375:n2233 | Required since 2022 by most epidemiology journals |
| Bonferroni for pheWAS-MR | Standard | Many outcomes; FDR if exploratory |
