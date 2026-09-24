# MRlap: sample overlap, winner's curse and weak instruments

Read when exposure and outcome GWAS may share samples (UKB-on-UKB, FinnGen-on-FinnGen), only summary statistics are available, or winner's curse is suspected. Verbatim from SKILL.md; the overlap decision rule is in "One-Sample vs Two-Sample Bias Direction". The winner's curse block at the end continues SKILL.md "Winner's curse at P~5e-8" (its heading, trigger and symptom stay there).

### MRlap: unified correction for sample overlap + winner's curse + weak instruments

MRlap (Mounier & Kutalik 2023 Genet Epidemiol 47:314) jointly corrects three biases that previously required three separate tools: sample overlap, winner's curse, and weak-instrument bias. It builds on an LDSC scaffold (cross-trait LD-score regression intercept estimates the overlap-induced covariance) and reweights the IVW estimate against the analytical bias-correction formula.

```r
remotes::install_github('n-mounier/MRlap')   # never on CRAN; bioconductor unsuitable
library(MRlap)

fit <- MRlap(
    exposure = gwas_X_df, exposure_name = 'BMI',
    outcome = gwas_Y_df, outcome_name = 'T2D',
    ld = 'eur_w_ld_chr/', hm3 = 'w_hm3.snplist',   # LDSC reference files
    MR_threshold = 5e-8, MR_pruning_dist = 500, MR_pruning_LD = 0.05
)
fit$MRcorrection$corrected_effect       # overlap + winner's curse + weak-IV corrected
fit$MRcorrection$corrected_effect_se
fit$LDSC$h2_exp                          # exposure heritability sanity check
fit$LDSC$int_crosstrait                  # cross-trait LDSC intercept; ~0 means no sample overlap
```

**Decision rule -- prefer MRlap when:** (a) any sample overlap is suspected, (b) only sumstats are available (no individual-level data for re-running GWAS on disjoint samples), (c) exposure discovery and outcome were both run inside the same biobank (UKB-on-UKB, FinnGen-on-FinnGen). MRlap returns NA / unstable estimates when h^2 < 0.05; in that regime, fall back to Burgess 2016 overlap-corrected IVW plus MR-RAPS for the weak-IV component.

### Winner's curse at P~5e-8 (Mechanism and Fix, moved from SKILL.md)

**Mechanism:** SNPs that just cross 5e-8 in discovery have over-estimated effect sizes (regression toward the mean in independent replication); MR uses inflated `beta_X`, biasing causal estimate.

**Fix:** (1) Three-sample design (discovery / replication-for-instrument-effect / outcome) where feasible. (2) When sumstats-only: MRlap (Mounier 2023), MR-SimSS (sample-splitting from sumstats), or RIVW (Ma 2023 Ann Statist 51:211 -- rerandomized IVW) jointly correct winner's curse + weak IVs + overlap. (3) Jiang 2023 IJE 52:1209 empirical magnitude: variant-level inflation ~50-400% near the genome-wide-significance threshold, dropping to <25% when the minimum P <= 1e-13.
