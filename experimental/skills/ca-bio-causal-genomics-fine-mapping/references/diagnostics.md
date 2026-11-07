# Critical LD Diagnostic Block (susie_rss)
**Goal:** Detect LD reference mismatch before reporting credible sets.

**Approach:** `estimate_s_rss()` quantifies the global Z-score / LD inconsistency as a scalar; `kriging_rss()` identifies individual SNPs whose Z-scores are inconsistent with the LD reference (typically genotyping errors, strand flips, or wrong reference panel).

```r
library(susieR)
s_hat <- estimate_s_rss(z = z_scores, R = ld_matrix, n = N)
# s_hat is the inferred scale of LD inconsistency.
# Source: susieR vignette "Diagnostic for summary statistic"; Zou 2022 PLoS Genet.
# Rule of thumb: s_hat < 0.05 acceptable; 0.05-0.10 marginal; > 0.10 refit or change LD reference.

cond_z <- kriging_rss(z = z_scores, R = ld_matrix, n = N)
# cond_z$conditional_dist is a data.frame (z, condmean, condvar, z_std_diff, logLR); flag abs(z_std_diff) > 3
# Common cause: strand flip, allele coding mismatch, or single-SNP imputation error.

# If diagnostic fails: refit with explicit scale parameter to absorb LD mismatch
fit <- susie_rss(z = z_scores, R = ld_matrix, n = N, L = 10, estimate_residual_variance = TRUE)
```

Skipping this block is the dominant cause of irreproducible fine-mapping. Always run before reporting credible sets.

## Per-Tool Failure Modes

### LD reference mismatch (most common)

**Trigger:** External LD matrix from 1000 Genomes / UK Biobank reference used for a GWAS conducted on a different cohort or ancestry mix.

**Mechanism:** Z-scores reflect the GWAS sample's LD; the reference R does not. The susie_rss likelihood depends on `z' R^{-1} z` being consistent with the modeled effects, and inconsistency manifests as spurious credible sets containing tag SNPs from the reference but not from the discovery cohort.

**Symptom:** `estimate_s_rss()` lambda > 0.05; `kriging_rss()` flags many SNPs with `|z_obs - z_exp| > 3`; credible sets contain physically distant SNPs (anti-correlated in LD with the lead) or include all SNPs at the locus.

**Fix:** Use in-sample LD whenever the cohort genotypes are accessible (compute with `plink --r square spaces` on the GWAS samples themselves — signed r, NOT `--r2`; the `spaces` modifier is required by FINEMAP, see `references/finemap-cli.md`). When only summary statistics are available, ancestry-stratify the LD reference exactly (e.g., 1000G EUR FIN+CEU+GBR+IBS+TSI for a Northern European GWAS, not full EUR). For mixed-ancestry GWAS, fine-map per ancestry then meta-analyze, or move to SuSiEx.

### Non-sparse architecture (biobank scale)

**Trigger:** Locus with one strong signal plus hundreds of weakly associated SNPs (polygenic shoulder); typical at biobank scale.

**Mechanism:** Vanilla SuSiE assumes a sparse sum-of-single-effects prior. With polygenic background, the model misallocates effects, producing inflated credible sets or many small spurious ones. Cui 2024 (Nat Genet 56:162) showed PIPs from SuSiE in this regime are systematically miscalibrated.

**Symptom:** Many small credible sets (5-15 per locus); replication in independent cohorts fails for non-lead credible sets; PIP distribution has a heavy tail.

**Fix:** Use SuSiE-inf or FINEMAP-inf (Cui 2024). These augment the sum-of-single-effects with an infinitesimal random-effect component that absorbs polygenic background. Source: github.com/FinucaneLab/fine-mapping-inf.

### L too small

**Trigger:** Locus with > 5 independent signals (HLA region, APOC1/APOE, LPA, IL6R region for some traits).

**Mechanism:** SuSiE assumes at most L independent effects. When the true number exceeds L, some signals are absorbed into existing components, distorting PIPs and credible sets for the captured signals.

**Symptom:** `length(fit$sets$cs)` equals L (all L slots used); credible set purity for higher-indexed sets is low (`fit$sets$purity[,'min.abs.corr'] < 0.5`); fits with larger L change top-PIP variants.

**Fix:** Increase L iteratively (L=10 -> 20 -> 30) until `length(fit$sets$cs)` < L (susieR auto-prunes unsupported effects so the returned CS count is the effective L). For HLA, start at L=30. The cost is mostly computational, not statistical: SuSiE prunes unused slots, so L=30 is safe when L=10 was right.

### prior_weights vs prior_variance confusion (PolyFun integration)

**Trigger:** Passing PolyFun output to susie_rss with `prior_variance=polyfun_priors` (wrong argument).

**Mechanism:** `prior_variance` in susie_rss is a single scalar (or vector of length L) for the per-effect variance, NOT a per-SNP probability. `prior_weights` is the per-SNP causal probability vector (sums to ~1). Passing PolyFun's per-SNP prior to `prior_variance` is silently accepted but applies a numerically nonsensical per-effect variance.

**Symptom:** PIPs nearly identical to the uniform-prior fit; functional annotations appear to have no effect.

**Fix:** Use `prior_weights = polyfun_priors$SNPVAR` (the PolyFun output column is uppercase `SNPVAR`; R is case-sensitive). Verify with `?susie_rss` in the installed version. Reference: github.com/omerwe/polyfun README, Weissbrod 2020 supplementary methods.

### Credible-set misinterpretation

**Trigger:** Reporting per-variant PIP without distinguishing "in credible set" from "high PIP".

**Mechanism:** The 95 percent credible-set guarantee is `P(causal variant in set) >= 0.95`. Per-variant PIPs within a set do not necessarily sum to 1 across all variants, and PIPs across overlapping sets can double-count posterior mass.

**Symptom:** Reporting "the top PIP variant" when the credible set is wide (size > 50); claiming a single variant is causal when the set contains 30 high-LD SNPs.

**Fix:** Always report (a) number of credible sets, (b) size of each set, (c) purity (`fit$sets$purity[,'min.abs.corr']`), (d) the top PIP variant within the set as the candidate lead. The credible set is the unit of inference; the top PIP variant is a candidate, not a conclusion.

### Cross-ancestry with single-ancestry LD

**Trigger:** Multi-ancestry meta-analyzed GWAS, then susie_rss with EUR LD.

**Mechanism:** Meta-analysis z-scores reflect a weighted mix of population LD structures; no single-population LD matrix matches.

**Fix:** Move to SuSiEx (joint cross-ancestry SuSiE; Yuan 2024). Per-ancestry fine-mapping followed by manual merging loses the shared-causal-variant information that SuSiEx exploits.

### Case-control GWAS passing Ntotal instead of Neff

**Trigger:** Passing `n = N_total` to `susie_rss()` for case-control GWAS derived from logistic regression.

**Mechanism:** susie_rss expects the effective sample size that determined the standard errors. For case-control logistic regression, `Neff = 4 / (1/Ncase + 1/Ncontrol)`; when cases are rare, total N can exceed Neff by 25x or more. Passing Ntotal rescales z-scores into a regime SuSiE never sees and makes the implied prior variance wrong.

**Symptom:** PIPs systematically biased; credible sets either too narrow (PIPs collapse to a single SNP that is not robust) or too wide (PIPs flatten); replication poor; sometimes z-score scale warnings from susieR.

**Fix:** `Neff = 4 / (1/Ncase + 1/Ncontrol)`. Example: Ncase=5000, Ncontrol=495000 -> Neff ~= 19,800 (NOT 500,000). For quantitative traits from linear regression, `n = N_total` is correct. Reference: Privé F et al 2022 HGG Adv 3:100136 (`bigsnpr` documents Neff handling); Willer 2010 Bioinformatics (METAL Neff convention).
