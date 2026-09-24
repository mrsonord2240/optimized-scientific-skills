---
name: bio-causal-genomics-fine-mapping
category: Data Analysis
description: Resolves GWAS associations to candidate causal variants and credible sets via SuSiE, susie_rss, FINEMAP, CAVIAR, DAP-G, PAINTOR, PolyFun, SuSiEx, MultiSuSiE, and FOCUS. Use when narrowing a GWAS lead SNP to a 95 percent credible set, choosing between in-sample and reference LD, calibrating non-sparse loci with SuSiE-inf or FINEMAP-inf, integrating functional priors via PolyFun, fine-mapping across ancestries with SuSiEx, diagnosing LD mismatch via estimate_s_rss and kriging_rss, handling HLA or long-range LD, or feeding credible sets into coloc.susie for colocalization.
tool_type: r
primary_tool: susieR
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: susieR 0.14.2, coloc 5.2.3, FINEMAP 1.4.2 (bioconda), PolyFun (head of `omerwe/polyfun` 2024), PAINTOR V3.0, SuSiEx 1.1.2 (bioconda), DAP-G (`xqwen/dap` @ `875ba40`), pyfocus 0.8+, R 4.4.3, PLINK 1.9.0-b.8 / 2.0 alpha 7.6.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('susieR')` then `?susie_rss` to confirm argument names (e.g., `prior_weights` vs `prior_variance` semantics)
- CLI: `finemap --help`, `SuSiEx --help`, `PAINTOR --help`, `dap-g --help` to confirm flags
- Python: `polyfun.py --help`

If a call throws an error about an argument that no longer exists, introspect the installed function and adapt rather than retrying.

# Fine-Mapping

**"Narrow my GWAS locus to the variants likely to be causal"** -> Fit a sparse Bayesian regression that propagates LD into posterior inclusion probabilities (PIPs) and credible sets, then validate that credible sets correspond to physically reasonable haplotypes given the LD reference.

- R (summary statistics + LD): `susieR::susie_rss(z, R, n, L=10)` + `estimate_s_rss` LD diagnostic
- R (individual-level genotypes): `susieR::susie(X, y, L=10)`
- CLI (shotgun stochastic search): `finemap --sss --in-files master.z --n-causal-snps 5 --prob-conv-sss-tol 0.001` (checked on FINEMAP 1.4.2 — flag names differ from older docs; see `references/finemap-cli.md`)
- CLI (cross-ancestry joint): `SuSiEx --sst_file=eur.sst,eas.sst --n_gwas=N1,N2 --ref_file=eur.bim,eas.bim --ld_file=eur_ld,eas_ld --chr_col=1,1 --snp_col=2,2 --bp_col=3,3 --a1_col=4,4 --a2_col=5,5 --eff_col=6,6 --se_col=7,7 --pval_col=8,8 --chr=<chr> --bp=<start,end> --out_dir=<dir> --out_name=<name> --plink=<path/to/plink>` (checked on SuSiEx 1.1.2 — `--plink` is required and `--ld_file` is an output prefix, not a pre-built matrix; see `references/susiex-cross-ancestry.md`)
- Python (functional priors): `polyfun.py --compute-h2-L2` -> per-SNP priors -> susie_rss with `prior_weights=`
- Python (TWAS fine-mapping): `focus finemap` on gene-level Z-scores

**Platform note:** FINEMAP, SuSiEx, PAINTOR, and DAP-G are POSIX (Linux/macOS) command-line binaries with no native Windows build. On Windows, run them under WSL, or use `susie_rss` / SuSiE-inf as the equivalent inference path.

Fine-mapping is a Bayesian model selection problem; LD is not noise but structured prior information. Most failure modes trace back to one of three issues: (a) LD reference mismatched to the GWAS sample; (b) the sparse-effects prior being wrong for the locus (polygenic background); or (c) too small an L cap. The `estimate_s_rss()` lambda and `kriging_rss()` per-SNP diagnostic catch (a) before downstream credible sets are reported.

## Tool Install Notes

```r
install.packages('susieR')                # CRAN; pin >= 0.12.27 for stable susie_rss API
install.packages('coloc')                 # CRAN; >= 5.2.3 for coloc.susie
install.packages(c('ggplot2', 'patchwork', 'dplyr', 'readr'))
```

```bash
# FINEMAP (CLI binary, not an R package). christianbenner.com/finemap.me have no working
# static download link on this platform -- bioconda mirrors the real binary (verified: 1.4.2):
conda install -c bioconda finemap

# PolyFun (Python)
git clone https://github.com/omerwe/polyfun
pip install -r polyfun/requirements.txt
# Pre-baked baseline-LF priors: https://data.broadinstitute.org/alkesgroup/UKBB_LD/baselineLF2.2.UKB.tar.gz
# polyfun.py's pd.read_csv(..., delim_whitespace=True) calls need pandas<3 (kwarg removed in 3.0).

# PAINTOR (C++ CLI; needs NLopt + Eigen)
git clone https://github.com/bogdanlab/PAINTOR_V3.0
make

# SuSiEx (C++ CLI); bioconda has a working prebuilt binary (verified: 1.1.2):
conda install -c bioconda susiex
# or from source: git clone https://github.com/getian107/SuSiEx && make -C src

# DAP-G (CLI; needs GSL + OpenMP)
git clone https://github.com/xqwen/dap && cd dap/dap_src && make    # verified: builds clean

# FOCUS (Python; for TWAS fine-mapping)
pip install pyfocus

# PLINK 1.9 / 2.0 for LD matrix generation
conda install -c bioconda plink plink2
```

## Decision Tree by Experimental Scenario

| Scenario | Recommended workflow | Why |
|----------|---------------------|-----|
| Individual-level genotypes available (UKB, in-house cohort) | `susie(X, y, L=10)` | In-sample LD is exact; no mismatch fragility |
| Summary statistics only, ancestry matches reference panel | `susie_rss(z, R, n, L=10)` + `estimate_s_rss` diagnostic; harmonize alleles first (`references/allele-harmonization.md`) | Standard external-LD pattern; verify lambda < 0.05 |
| Single-locus EUR GWAS, sparse architecture | susie_rss with L=10, baseline functional priors optional; independent confirmation with FINEMAP (`references/finemap-cli.md`) | Most-common setting; SuSiE default works |
| Locus with strong polygenic shoulder (biobank scale) | SuSiE-inf (Cui 2024) | Adds infinitesimal component; calibrates non-sparse PIPs |
| Multi-ancestry GWAS (EUR + EAS + AFR) | SuSiEx with per-pop sumstats and LD (`references/susiex-cross-ancestry.md`) | Joint inference shrinks credible sets; per-ancestry meta loses LD information |
| Locus with > 5 expected independent signals (HLA, lipid loci) | susie_rss with L=20-30 | Default L=10 caps signal count; HLA needs extension |
| TWAS hits with co-regulated genes | FOCUS / MA-FOCUS | Variant-level fine-mapping cannot distinguish co-regulated gene candidates |
| Want functional priors (coding, conserved, regulatory) | PolyFun -> susie_rss with `prior_weights` (`references/polyfun-functional-priors.md`) | Genome-wide SLDSC priors sharpen PIPs more than locus-level annotations |
| QTL fine-mapping (eQTL, sQTL, caQTL) at transcriptome scale | DAP-G + TORUS OR susie_rss per gene (DAP-G command: `references/dap-g-cli.md`; tool comparison: `references/method-comparison.md`) | DAP-G is built for QTL throughput; SuSiE works per gene |
| Low-N QTL (GTEx tissue panel, N < 1000) | susie_rss with `coverage = 0.9` (or 0.8); document choice | Default 0.95 returns very wide credible sets at low power; report the relaxed coverage explicitly in methods |
| HLA region (chr6:28-34 Mb) or chr8 inversion | Specialized workflow: stratify haplotypes; consider HLA-specific imputation; or exclude (`references/hla-long-range-ld.md`) | LD structure is too complex; standard methods unreliable |
| Cross-feed into colocalization | susie_rss -> coloc.susie() (`references/coloc-susie.md`) | Modern coloc operates on credible sets, not single SNPs |

## Reference Files

Read only the file the request needs.

| File | Read when |
|------|-----------|
| `references/method-comparison.md` | Choosing between SuSiE, FINEMAP, CAVIAR, DAP-G, PAINTOR, PolyFun, SuSiEx, MultiSuSiE, FOCUS (Algorithmic Taxonomy), or two methods disagree (Reconciliation table) |
| `references/allele-harmonization.md` | Sumstats and LD reference come from different sources, or `estimate_s_rss` / `kriging_rss` flag many SNPs with allele-swap symptoms |
| `references/polyfun-functional-priors.md` | Functional priors: PolyFun commands, `SNPVAR` -> `prior_weights`, manual coding-variant priors |
| `references/susiex-cross-ancestry.md` | Multi-ancestry GWAS; `SuSiEx` command and population-order rule |
| `references/finemap-cli.md` | Independent FINEMAP confirmation: master file and `--sss` run |
| `references/dap-g-cli.md` | QTL-scale fine-mapping or an independent DAP-G confirmation: `-d_z`/`-d_ld` command and output format |
| `references/coloc-susie.md` | Passing credible sets to `coloc.susie` (SNP-name precondition and guard) |
| `references/hla-long-range-ld.md` | HLA, chr8 inversion or other long-range LD locus |
| `references/reviewer-pushback.md` | Writing methods or a rebuttal: standard responses to fine-mapping reviewer questions |

## Critical LD Diagnostic Block (susie_rss)

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

## Reconciliation: When Methods Disagree

The table of disagreement patterns (likely cause and action) is in `references/method-comparison.md`.

**Operational rule:** For high-confidence reporting, require that (a) `estimate_s_rss()` lambda < 0.05; (b) at least one credible set has purity > 0.5 (`min_abs_corr >= 0.5`, equivalent to r2 >= 0.25); (c) the lead PIP variant within that set is reproduced by an independent method (FINEMAP, SuSiEx, or in-sample SuSiE if reference-LD was used). Anything failing these three is exploratory.

## Quantitative Thresholds

| Quantity | Threshold | Source / Rationale |
|----------|-----------|-------------------|
| Credible set coverage (well-powered GWAS) | 0.95 (default) | Wang 2020 JRSSB; standard convention |
| Credible set coverage (low-N eQTL, GTEx tissue) | 0.9 or 0.8 | At N < 1000, default 0.95 returns very wide CS; document choice in methods |
| Credible set purity (rare-variant fine-mapping) | min_abs_corr >= 0.1 | LD genuinely sparse; relax to retain signal |
| Credible set purity (default common-variant) | min_abs_corr >= 0.5 (r2 >= 0.25) | susieR default; below this the set is LD-confounded |
| Credible set purity (publication-strict) | min_abs_corr >= 0.7 | Stringent claim; rare in practice |
| PIP suggestive | > 0.5 | Convention; "more likely than not causal among set" |
| PIP strong | > 0.9 | Convention; high-confidence single candidate |
| PIP very strong | > 0.95 | Convention; near-certain candidate within credible set |
| L (default cap) | 10 | susieR default; sufficient for most non-HLA loci |
| L (HLA / complex loci) | 20-30 | Empirical; HLA hosts > 10 independent signals for many traits |
| `n` for case-control susie_rss | Neff = 4/(1/Ncase + 1/Ncontrol), NOT Ntotal | Privé F et al 2022 HGG Adv 3:100136; matches the SE scale of logistic-regression sumstats |
| `estimate_s_rss` lambda acceptable | < 0.05 | susieR vignette; > 0.10 indicates serious LD mismatch |
| `kriging_rss` per-SNP flag | abs(`z_std_diff`) > 3 (`|z_obs - z_exp|` elsewhere in this file) | susieR vignette; flag for manual review |
| Locus window (default) | +/- 500 kb from sentinel | Conventional; covers most LD blocks |
| Locus window (conditional-p floor) | Extend until conditional -log10(p) < 4 | Avoids truncating a secondary signal whose conditional evidence leaks into the window edge |
| Locus window (long-range LD) | 5+ Mb or stratify | HLA chr6:25-35Mb, chr8 inversion chr8:8.1-11.9Mb hg38, chr17 H1/H2 inversion |
| FINEMAP `--n-causal-snps` | 5 | Default; raise for HLA |
| FINEMAP `--prob-conv-sss-tol` | 0.001 | Convergence tolerance; rarely needs change (checked on FINEMAP 1.4.2) |

## TWAS Fine-Mapping (FOCUS) -- delegated

Variant-level fine-mapping cannot distinguish causal genes among co-regulated TWAS hits. FOCUS / MA-FOCUS extend fine-mapping to the predicted-expression level; see causal-genomics/transcriptome-wide-association for the full FOCUS workflow and reconciliation with variant-level credible sets.

## Required Reporting Schema for Fine-Mapping

Every locus reported should carry these columns; missing fields are the most common reviewer complaint.

| Column | Description |
|--------|-------------|
| locus_id | Locus identifier (chr:start-end or sentinel rsID) |
| method | susie_rss / FINEMAP / PAINTOR / SuSiEx / SuSiE-inf |
| L_used | `sum(fit$V > 0)` (effective L; not just the cap passed in; susieR sets pruned effects' prior variance to 0) |
| n_credible_sets | Number of returned credible sets at the chosen coverage |
| cs_size | Variants per credible set |
| cs_purity_min / cs_purity_mean | min and mean `fit$sets$purity[,'min.abs.corr']` per set |
| top_pip_snp | Lead variant in each credible set |
| top_pip | Posterior inclusion probability of top_pip_snp |
| lambda_s | `estimate_s_rss` diagnostic for the locus |
| kriging_outlier_count | Count of SNPs with `|z_obs - z_exp| > 3` |
| ld_panel | 1KG-EUR / UKB-EUR / in-sample / TopMed |
| prior_source | uniform / PolyFun-EUR / PolyFun matched-ancestry baseline-LF / manual coding-variant |
| coverage | 0.95 default; 0.9 or 0.8 documented for low-N |
| n_effective | Sample size passed to susie_rss (Neff for case-control) |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `IBSS algorithm did not converge` warning | L too small OR LD mismatch | Increase L; run `estimate_s_rss`; check ancestry match |
| Credible set contains all SNPs at locus | LD reference matches discovery poorly; `s_hat` > 0.1 | Switch to in-sample LD or stratify reference ancestry |
| PIPs identical to GWAS p-value rank | Effectively no LD information used; check LD matrix orientation | Verify SNP order in z and R match exactly; check for transposed R |
| Negative eigenvalues in LD matrix | Numerical PSD violation from finite-precision storage | Add small ridge: `R <- R + diag(1e-4, nrow(R))`; or use `Matrix::nearPD` |
| `pip` all ~ 1/p (uniform) | Convergence failure OR all effects pruned | Check `fit$converged`; raise L; check Z scale |
| FINEMAP `Error: SNP names do not match` | .z and .ld SNP order differ | Ensure both are sorted identically; pass matched .snp file |
| Coloc.susie returns NULL | One trait has zero credible sets | Verify both fits succeeded; lower coverage to 0.9 if signal is weak |
| Coloc.susie crashes with `data.table` error (`:= is defined for use in j`) | `z`/`R` passed to `susie_rss` without SNP-ID names, so `coloc.susie`'s internal SNP match is empty | Name `z1`/`z2` and set matching `ld_matrix` dimnames before fitting; see the precondition in `references/coloc-susie.md` |
| SuSiEx output empty, or fails inside its own internal PLINK calls (`Error: No variants remaining after --extract`) | Per-population lists misaligned with reference panels, or `--plink=<path>` omitted (required on SuSiEx 1.1.2, not caught at argument-parse time) | Verify `--sst_file`/`--ref_file`/`--ld_file` share the same population order and `--bp` window; pass `--plink=$(command -v plink)` — see `references/susiex-cross-ancestry.md` |
| FINEMAP `Cannot recognize flag` (`--prob-tol`/`--n-iterations`/`--n-convergence`); `.ld` "Expected N SNPs ... encountered only 1"; or no plain `locus.cred` written | Flags renamed in FINEMAP 1.4.2; `plink --r square` defaults to tab-delimited (FINEMAP needs spaces); output is `<prefix>.cred<k>` per causal-count model, not a plain `.cred` | Use `--prob-conv-sss-tol`/`--n-iter`/`--n-conv-sss`; add `--r square spaces`; read `locus.cred1`, `locus.cred2`, ...; see `references/finemap-cli.md` |
| `dap-g` exits with status 1 | Normal — DAP-G returns 1 even on a fully successful run | Judge success by stdout (`Independent association signal clusters`), never the exit code |
| PolyFun priors do not change PIPs | Passed to `prior_variance` instead of `prior_weights` | Read susieR docs; use `prior_weights=` |

## References

- Wang G, Sarkar A, Carbonetto P, Stephens M 2020 J R Stat Soc B 82:1273 (SuSiE / IBSS)
- Zou Y, Carbonetto P, Wang G, Stephens M 2022 PLoS Genet 18:e1010299 (susie_rss for summary statistics)
- Cui R, Elzur RA, Kanai M, Ulirsch JC, Weissbrod O et al 2024 Nat Genet 56:162 (SuSiE-inf / FINEMAP-inf for non-sparse loci)
- Benner C, Spencer CC, Havulinna AS, Salomaa V, Ripatti S, Pirinen M 2016 Bioinformatics 32:1493 (FINEMAP)
- Hormozdiari F, Kostem E, Kang EY, Pasaniuc B, Eskin E 2014 Genetics 198:497 (CAVIAR)
- Wen X, Lee Y, Luca F, Pique-Regi R 2016 AJHG 98:1114 (DAP-G)
- Kichaev G, Yang WY, Lindstrom S, Hormozdiari F, Eskin E et al 2014 PLoS Genet 10:e1004722 (PAINTOR)
- Weissbrod O, Hormozdiari F, Benner C, Cui R, Ulirsch J et al 2020 Nat Genet 52:1355 (PolyFun + functional priors)
- Yuan K, Longchamps RJ, Pardinas AF, Yu M, Chen TT et al 2024 Nat Genet 56:1841 (SuSiEx cross-ancestry)
- Rossen J, Shi H, Strober BJ, Zhang MJ, Kanai M, McCaw ZR, Liang L, Weissbrod O, Price AL 2025 Nat Genet 58:67 (MultiSuSiE; doi:10.1038/s41588-025-02450-5)
- Mancuso N, Freund MK, Johnson R, Shi H, Kichaev G et al 2019 Nat Genet 51:675 (FOCUS for TWAS fine-mapping)
- Wallace C 2021 PLoS Genet 17:e1009440 (coloc.susie integration)
- Schaid DJ, Chen W, Larson NB 2018 Nat Rev Genet 19:491 (fine-mapping review)
- Hutchinson A, Asimit J, Wallace C 2020 Hum Mol Genet 29:R81 (fine-mapping review)

## Related Skills

- causal-genomics/colocalization-analysis - coloc.susie consumes susie_rss credible sets; equivalent harmonize helper
- causal-genomics/effector-gene-prioritization - Downstream gene-assignment from credible-set variants
- causal-genomics/transcriptome-wide-association - FOCUS / MA-FOCUS for gene-level fine-mapping
- causal-genomics/genomic-sem - Joint multi-trait fine-mapping when credible sets are shared across traits
- causal-genomics/mendelian-randomization - Fine-mapped variants as cis-instruments
- causal-genomics/pleiotropy-detection - Per-credible-set pleiotropy testing
- atac-seq/enhancer-gene-linking - ABC / ENCODE-rE2G linking credible-set variants to target genes
- population-genetics/linkage-disequilibrium - Constructing LD matrices for susie_rss
- population-genetics/association-testing - Upstream GWAS summary statistic generation
- workflows/gwas-pipeline - End-to-end GWAS pipeline producing fine-mapping input
- variant-calling/variant-annotation - Annotating credible-set variants with VEP / coding consequence
- pathway-analysis/go-enrichment - Downstream gene-level interpretation of credible-set targets
