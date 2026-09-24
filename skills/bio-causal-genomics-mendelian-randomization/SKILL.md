---
name: bio-causal-genomics-mendelian-randomization
category: Data Analysis
description: Estimate causal effects of an exposure on an outcome from GWAS summary statistics using genetic instruments. Implements IVW (fixed/random), MR-Egger, weighted median/mode, MR-RAPS, CAUSE, GSMR-HEIDI, MR-PRESSO, MVMR, MR-Clust, LCV, and LHC-MR via TwoSampleMR, MendelianRandomization, MR-PRESSO, cause, and lhcMR. Use when testing causal direction between traits, evaluating drug-target effects via cis-pQTL/cis-eQTL, performing multivariable mediation MR, distinguishing causation from correlated horizontal pleiotropy, or producing STROBE-MR-compliant sensitivity batteries.
tool_type: mixed
primary_tool: TwoSampleMR
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: TwoSampleMR 0.6.0+, MendelianRandomization 0.10+, MR-PRESSO 1.0+, cause 1.2+, MVMR 0.4+, ieugwasr 1.0+, MRlap 0.0.3.2+, coloc 5.2+, mrclust 0.1+, lhcMR 0.0.1+, R 4.4+. Both TwoSampleMR 0.6.0 and ieugwasr 1.0 are the JWT-transition versions; older versions still expect deprecated OAuth.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI (plink, GCTA-GSMR): `<tool> --version` then `<tool> --help`

If code throws an error referencing a function that has moved (e.g. `ieugwasr::ld_clump` vs `TwoSampleMR::clump_data`) or an OAuth token failure, introspect the installed API and adapt the example rather than retrying.

# Mendelian Randomization

**"Test whether trait X causally affects trait Y from GWAS summary statistics"** -> Use genetic variants as instrumental variables (IVs) that satisfy three assumptions (relevance, independence, exclusion restriction) to estimate `beta_causal = beta_outcome / beta_exposure` under the IV framework (Davey Smith & Ebrahim 2003 IJE 32:1; Burgess & Thompson 2021 Chapman & Hall/CRC, 2nd ed.). Tool choice is a decision about the **regime** (one-sample vs two-sample, sparse vs polygenic, drug-target vs polygenic exposure) and the **pleiotropy model** (balanced, directional InSIDE, correlated horizontal). Wrong tool inflates Type-I error or attenuates true effects in a direction predictable from the bias structure.

- R: `TwoSampleMR::mr()` orchestrates IVW + Egger + weighted median + weighted mode in one call
- R: `MendelianRandomization::mr_ivw / mr_egger / mr_median / mr_mbe / mr_conmix` per-method API (S4 objects; MR-RAPS is NOT in this package -- use `TwoSampleMR::mr_raps()` which wraps the GitHub `mr.raps`)
- R: `MRPRESSO::mr_presso()` global / outlier / distortion tests
- R: `cause::cause()` correlated horizontal pleiotropy mixture
- R: `MVMR::strength_mvmr() + MVMR::ivw_mvmr()` multivariable conditional-F + IVW

## Reference Files

Read a file only when the request needs that method; everything every request needs (scope, decision tree, thresholds, standard workflow, Common Errors, install) stays in this file.

| File | Read when |
|------|-----------|
| `references/mr-presso.md` | Running MR-PRESSO, extracting outlier SNPs, choosing `NbDistribution` |
| `references/mvmr-conditional-f.md` | Any multivariable MR or mediation MR; conditional F, `qhet_mvmr` fallback |
| `references/mrlap-overlap-correction.md` | Sample overlap suspected (UKB-on-UKB), or sumstats-only correction of winner's curse + weak IV (includes the winner's curse mechanism and fix) |
| `references/simex-egger-nome.md` | `I^2_GX < 0.9` and Egger must be reported |
| `references/bidirectional-steiger.md` | Reverse MR, Steiger filtering or directionality, a Steiger flag that contradicts biology |
| `references/cis-mr-and-binary-outcomes.md` | Drug-target cis-MR, binary or case-only outcomes, collider bias in stratified MR |
| `references/cohorts-and-reviewer-pushback.md` | Choosing a biobank/ancestry, or answering a reviewer's overlap/weak-IV/pleiotropy/reverse-causation/InSIDE critique |
| `references/tool-installation.md` | A package is missing or fails to install; setting up local clumping |
| `references/bibliography.md` | Citing a method |

Runnable code lives in `scripts/` (run with `Rscript`): `twosample_workflow.R` (standard workflow, writes `harmonised.tsv`), `mr_presso_outliers.R`, `mvmr_conditional_f.R`, `simex_egger.R`. Each has a usage header.

## Statistical Model Taxonomy

| Method | Pleiotropy assumption | Min instruments | Strength | Fails when |
|--------|------------------------|-----------------|----------|------------|
| IVW (fixed) | All IVs valid | 2 | Most efficient under no pleiotropy | Any directional or balanced pleiotropy inflates Type-I |
| IVW (random effects) | Balanced + InSIDE | 3 | Standard primary; absorbs heterogeneity into wider SE | Directional pleiotropy biases the point estimate |
| MR-Egger | Directional pleiotropy + InSIDE | 10+ for power | Detects + corrects directional pleiotropy via intercept (Bowden 2015 IJE 44:512) | NOME violated (`I^2_GX < 0.9`); SIMEX correction required; underpowered <10 SNPs |
| Weighted median | Up to 50% invalid IVs | 3 | Robust to a minority of bad instruments (Bowden 2016 Genet Epidemiol 40:304) | >50% invalid IVs |
| Weighted mode | Zero modal pleiotropy (ZEMPA) | 3 | Robust if the modal estimate is unbiased (Hartwig 2017 IJE 46:1985) | Bimodal pleiotropy; small numbers |
| MR-RAPS | Balanced pleiotropy + weak instruments | 10+ | Profile-score robust to weak-IV + balanced horizontal pleiotropy (Zhao 2020 Ann Stat 48:1742) | Strong directional pleiotropy; CRAN-archived 2025-03-01 |
| CAUSE | Correlated horizontal pleiotropy (CHP) | 100+ sig SNPs | Explicit shared-factor mixture; protects against CHP-driven false positives (Morrison 2020 Nat Genet 52:740) | Sparse polygenic exposures; <100 sig SNPs |
| GSMR + HEIDI-outlier | Outlier removal under InSIDE | 10+ | Alternative outlier detection; integrates with LD reference (Zhu 2018 Nat Commun 9:224) | Requires individual-level LD; HEIDI conservative |
| MR-PRESSO | Outlier-driven horizontal pleiotropy | 4+ | Global / outlier / distortion three-step (Verbanck 2018 Nat Genet 50:693) | Blind to CHP; computationally heavy at large NbDistribution |
| MVMR (IVW) | Conditional independence after measured pleiotropy | 1+ per exposure | Accounts for measured horizontal pleiotropy via multivariable regression (Sanderson 2019 IJE 48:713) | Conditional F < 10 on any exposure |
| MR-Clust | Heterogeneous causal effects (multiple mechanisms) | 30+ | Clusters SNPs by their causal-effect estimate (Foley 2021 Bioinformatics 37:531) | Single causal mechanism; small instrument sets |
| Contamination mixture | Mixture of valid + invalid IVs | 10+ | Profile-likelihood mixture (Burgess 2020 Nat Commun 11:376) | Sparse signal |
| LCV | Genome-wide; distinguishes causation vs genetic correlation | All SNPs | Tests `gcp` parameter using LDSC-style block jackknife (O'Connor & Price 2018 Nat Genet 50:1728) | Two-trait covariance dominated by a third confounder |
| LHC-MR | Bidirectional + heritable confounder | All SNPs | Joint likelihood over genome-wide markers; estimates both directions + confounder (Darrous 2021 Nat Commun 12:7274) | Computationally heavy; rare-variant trait |
| MRlap | Sample overlap + winner's curse + weak-IV jointly | Genome-wide sumstats | LDSC-scaffolded joint correction (Mounier & Kutalik 2023 Genet Epidemiol 47:314) | LDSC intercept poorly estimated (h^2 < 0.05); non-EUR without matched LD scores |
| Doubly-Ranked MR (DRMR) | Non-linear, non-parametric | 5+ strata | Non-parametric stratification (Tian 2023 PLoS Genet 19:e1010823); replaces residual stratification when linearity fails | Continuous exposures only; needs individual-level data; Hamilton 2023 medRxiv 23293658 shows stratum-specific bias from age/sex |

Methodology evolves; benchmark consensus shifts every 2-3 years. Verify against the current Slob & Burgess 2020 *Genet Epidemiol*, Burgess 2023 *Wellcome Open Res* "Guidelines for performing Mendelian randomization" (v3+), and STROBE-MR 2021 reporting standards before locking a method as primary.

## Decision Tree by Experimental Scenario

| Scenario | Primary method | Sensitivity battery | Why |
|----------|----------------|----------------------|-----|
| Standard two-sample, independent cohorts, polygenic exposure | IVW (random) | Egger + weighted median + MR-PRESSO + MR-RAPS + Steiger | Default; covers balanced, directional, outlier, weak-IV regimes; MR-PRESSO details in `references/mr-presso.md` |
| One-sample (e.g. UK Biobank both ends) | IVW, with MR-RAPS as sensitivity if mean F is borderline | Egger + LCV + jackknife SE + MRlap/Burgess-2016 overlap correction | One-sample F-stat floor shifts to F >= 20; jackknife SE preferred over analytic at one-sample scale; do NOT run exposure GWAS and outcome GWAS on the same individuals then claim two-sample (Barry 2021 PLoS Genet 17:e1009703 collider bias); within-stratum MR (e.g. "MR among smokers") risks collider bias from the stratification variable; MR-RAPS corrects weak-IV bias only and does NOT correct sample-overlap/confounding bias even when F is already high (see Operational rule below); overlap correction in `references/mrlap-overlap-correction.md` |
| Partial sample overlap (UKB exposure + UKB outcome) | MR-RAPS with overlap correction | Sample-overlap-adjusted IVW (Burgess 2016 Genet Epidemiol 40:597) | Bias is intermediate, proportional to z-score correlation; MRlap in `references/mrlap-overlap-correction.md` |
| Drug-target / cis-MR (cis-pQTL, cis-eQTL) | IVW restricted to cis window | Colocalization PP.H4 + LD-prune within window | Exclusion restriction relaxed because the protein/transcript directly mediates effect (Schmidt 2020 Nat Commun 11:3255); see `references/cis-mr-and-binary-outcomes.md` |
| MVMR for measured pleiotropy (e.g. LDL adjusted for HDL/TG) | `MVMR::ivw_mvmr` | Conditional F + Q_A heterogeneity | Required when exposures correlate via shared SNPs; code in `references/mvmr-conditional-f.md` |
| Mediation MR (X -> M -> Y) | MVMR difference of total vs direct | Two-step MR + product-of-coefficients (Carter 2021 Eur J Epidemiol 36:465) | Network MR; quantifies indirect effect; MVMR code in `references/mvmr-conditional-f.md` |
| Polygenic exposure with potential CHP (e.g. BMI -> CHD) | CAUSE (primary) + IVW (secondary) | Egger + MR-PRESSO + LCV | CAUSE explicitly models CHP via shared-factor; needs >=100 sig SNPs |
| Binary outcome (e.g. T2D) on linear scale | IVW on log-OR with log-additive coding | All sensitivity on log-OR; report exp(beta) | Linearity of MR estimating equation holds on log-OR not OR; non-collapsibility in `references/cis-mr-and-binary-outcomes.md` |
| Time-to-event (Cox) outcome | IVW on log-HR | Non-collapsibility-aware log-HR reporting | Non-collapsibility caveats apply |
| Non-linear MR (e.g. alcohol J-curve) | DRMR (Tian 2023) + residual stratification side-by-side | Negative-control outcomes (genotype-vs-sex within strata); Hamilton 2023 limitation cited | Both methods produce stratum-specific bias from age/sex effects (Hamilton 2023 medRxiv 23293658); pre-specify the non-linear hypothesis, do not data-snoop the J-curve, report negative-control sanity checks |
| Single-patient rare disease | Not MR -- use FRASER/DROP outlier framework | See alternative-splicing/outlier-splicing-detection | MR requires summary stats; n=1 is wrong regime |

## One-Sample vs Two-Sample Bias Direction

| Design | Weak-IV bias direction | Reason |
|--------|------------------------|--------|
| One-sample, F<10 | Toward confounded observational estimate (overestimates causal effect if confounding is in same direction) | Sample correlation between IV-X and IV-Y residuals |
| Two-sample non-overlapping, F<10 | Toward null | Independent samples decouple residuals (Burgess 2011 IJE 40:755) |
| Two-sample with partial overlap | Intermediate; proportional to overlap fraction and z-score correlation | Burgess 2016 Genet Epidemiol 40:597; correction available |

**Operational rule:** Whenever both GWAS came from UK Biobank (or any single biobank), treat the analysis as one-sample-equivalent. MR-RAPS corrects weak-instrument bias in this regime; it does NOT correct the separate sample-overlap/confounding bias that one-sample-equivalent designs also carry, and the two can coexist even when F is well above the one-sample floor (F >= 20) -- a real UKB-on-UKB run at mean F=44.6 with a planted null effect showed naive IVW biased (b=0.18) and MR-RAPS no better (b=0.196), because the bias source was sample-overlap confounding, not weak instruments. Use MR-RAPS when F is borderline; use MRlap (below) or Burgess 2016 overlap-corrected IVW whenever overlap/confounding is suspected, regardless of F. Treating it as "two-sample because separate GWAS files" is a common error and produces overestimates.

MRlap (joint overlap + winner's curse + weak-IV correction): code, decision rule and the h^2 < 0.05 fallback are in `references/mrlap-overlap-correction.md`. Prefer it whenever any sample overlap is suspected, and always for UKB-on-UKB.

## Drug-Target / cis-MR Framework

cis-MR restricts instruments to the cis window (+/-500 kb, clump r2 < 0.1, coloc PP.H4 >= 0.7, flag protein-altering variants); the full framework, plus binary outcomes, non-collapsibility and collider bias in case-only/stratified MR, is in `references/cis-mr-and-binary-outcomes.md`. Drug-target nominations belong to causal-genomics/proteome-mr-drug-target.

## Per-Method Failure Modes

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

## TwoSampleMR Standard Workflow

**Goal:** Produce a defensible primary IVW estimate plus a full sensitivity battery from two-sample summary statistics.

**Approach:** Extract genome-wide significant instruments -> clump (local plink preferred) -> extract outcome -> harmonise -> mr -> pleiotropy + heterogeneity + leave-one-out -> Steiger -> MR-PRESSO -> report.

```bash
Rscript scripts/twosample_workflow.R --exposure exposure_gwas.tsv --outcome outcome_gwas.tsv \
    --outdir mr_out --bfile 1kg_EUR/EUR    # --no-clump ONLY for instruments already LD-independent
```

`scripts/twosample_workflow.R` (TwoSampleMR, ieugwasr) reads two TSVs with columns `SNP, BETA, SE, A1, A2, EAF, P, N` (`N` is required: `directionality_test()` needs `samplesize_col`) and runs: `P < 5e-8` -> per-instrument F from the EXPOSURE `(beta/se)^2 >= 10` (Staiger-Stock 1997) -> `ld_clump()` at r2 = 0.001 / 10 Mb (local plink; `genetics.binaRies::get_plink_binary()` unless `--plink`) -> `harmonise_data(action = 2)` (infer from EAF; drops MAF~0.5 palindromes) -> `mr()` IVW / Egger / weighted median / weighted mode -> `mr_heterogeneity()` (Cochran Q), `mr_pleiotropy_test()` (Egger intercept), `mr_leaveoneout()` -> `directionality_test()`, which the script turns into an error when it silently returns NULL (missing p-value/sample-size columns; supply `r.exposure`/`r.outcome`, e.g. via `get_r_from_lor()`, for binary traits). It writes `primary`, `heterogeneity`, `pleiotropy`, `leaveoneout`, `steiger` and `harmonised` TSVs to `--outdir`; `harmonised.tsv` feeds the MR-PRESSO and SIMEX scripts.

## MR-PRESSO Outlier Detection

Global / outlier / distortion tests, the corrected outlier rule (adjusted P <= SignifThreshold, never divided by `nrow(dat)` again), the mr_keep filter and the `NbDistribution` cost/floor are in `references/mr-presso.md`. Read it before extracting outlier SNPs.

## CAUSE for Correlated Horizontal Pleiotropy

CAUSE (Morrison 2020 Nat Genet 52:740) fits a shared-factor mixture to genome-wide sumstats and compares causal vs sharing-only models by delta-ELPD. Workflow: `gwas_merge()` -> sample ~1M variants for `est_cause_params()` -> filter sig SNPs (P < 1e-3) and optionally LD-prune -> `cause(X, variants, param_ests)`. Needs >= 100 sig SNPs. Full annotated example and ELPD interpretation in causal-genomics/pleiotropy-detection. Version note: `cause` 1.2.0.335 with `loo` 2.10.1 crashes inside `cause:::in_sample_elpd_loo` ("non-numeric argument to binary operator", from a `loo_compare()` return-shape change), even with 147 significant SNPs; that is a package mismatch, not a data problem, so check `packageVersion('loo')` before debugging the input.

## MVMR with Conditional F

Format, conditional F per exposure (> 10 required; guard at F < 1), `ivw_mvmr`, Q_A and the `qhet_mvmr` fallback caveats are in `references/mvmr-conditional-f.md`. Read it for any multivariable MR.

## Bidirectional and Steiger

Reverse MR with an independent instrument set, `steiger_filtering()` (per SNP) and `directionality_test()` (global; returns NULL silently without sample sizes) with the operational reporting rule: `references/bidirectional-steiger.md`.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| IVW sig, Egger null with non-zero intercept | Directional pleiotropy | Report Egger as primary if `I^2_GX >= 0.9`; otherwise SIMEX-correct |
| IVW sig, weighted median sig, mode null | Mode underpowered or bimodal pleiotropy -- only if Cochran Q and the MR-PRESSO global test are non-significant | Trust the agreement of IVW + median only then; if Q or the PRESSO global test is significant, treat as the next row |
| IVW, weighted median and MR-RAPS all sig; Egger and mode null; Q and PRESSO global sig | Directional pleiotropy the median and RAPS do not absorb (planted null effect of 0, 30% invalid IVs: IVW p=2e-6, RAPS p=5e-5, median p=0.006, Egger p=0.72, mode p=0.82; outlier-corrected IVW still p=1.5e-4) | Do not report a causal effect; report as consistent with pleiotropy, lead with Egger/mode and the heterogeneity result |
| MR-PRESSO distortion test sig | Outliers materially shift estimate | Report PRESSO-adjusted estimate as primary |
| CAUSE sig, IVW sig, same direction | High-confidence causal claim | Report both; emphasize CAUSE rules out CHP |
| CAUSE null, IVW sig, same direction | CHP indistinguishable from causation | Downgrade to "consistent with causation but not separable from CHP" |
| Forward MR sig, reverse MR also sig | Bidirectional causation OR confounded | Use LHC-MR or LCV for genome-wide resolution |
| IVW sig but `mean F << 10` in one-sample design | Weak-IV bias toward observational | Re-run with MR-RAPS; report adjusted estimate |

**Operational rule for publication:** Primary IVW + concordant Egger (or weighted median if NOME violated) + non-significant MR-PRESSO global test + Steiger correct direction = publication-ready evidence. CAUSE concordance is required when the exposure is polygenic and the prior on CHP is high (e.g. BMI -> outcome, lipids -> outcome). Single-method "significant IVW" claims should be downgraded to exploratory.

## Cohort Gotchas and Anticipated Reviewer Pushback

Ancestry mixing in UKB/FinnGen/MVP/GBMI/AoU, UKB-on-UKB one-sample-equivalent bias, and the standard responses to the six usual reviewer pushbacks (overlap, weak instruments, pleiotropy, reverse causation, pre-registration, InSIDE): `references/cohorts-and-reviewer-pushback.md`.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| F-statistic computed from outcome | Reading off `beta.outcome / se.outcome` | Compute from EXPOSURE; outcome F is meaningless for IV strength |
| All instruments dropped at harmonise | All SNPs palindromic; or allele coding mismatch | Verify A1/A2 conventions match; try `action = 3` to drop, not assume |
| `Unauthorized` / `403` from OpenGWAS | OAuth deprecated May 2024; JWT token required | Generate at api.opengwas.io -> set `OPENGWAS_JWT=<token>` in `~/.Renviron` -> restart R -> verify with `ieugwasr::get_opengwas_jwt()`. For production, skip the API: use `ld_clump_local()` with local 1KG bfile (saves rate-limit + auth headaches). |
| Conditional F vs total F confusion in MVMR | Reporting `mean(F)` not `strength_mvmr()` per-exposure | Always use `MVMR::strength_mvmr()` and report each column |
| `install.packages("mr.raps")` fails | CRAN-archived 2025-03-01 | `remotes::install_github('qingyuanzhao/mr.raps')`; TwoSampleMR's `mr_raps()` wrapper internally calls this package |
| MR-PRESSO returns NA p-value, or warns "Outlier test unstable" | `NbDistribution` too small (needs > `nrow(dat) / SignifThreshold`); signal too thin | Increase to >= 5000 (>= 10000 for publication; run in the background); check that >= 4 SNPs remain after harmonization |
| MR-PRESSO flags no outliers though the global test is significant | Outlier threshold divided by `nrow(dat)` a second time (its outlier P is already Bonferroni-adjusted) | Use adjusted P <= `SignifThreshold` as in "MR-PRESSO Outlier Detection" |
| Egger intercept "highly significant" with 5 SNPs | Underpowered Egger over-fits the slope | Egger needs >= 10 SNPs; below that, intercept is unreliable |
| Sample-overlap correction ignored | Treating UKB-on-UKB as two-sample | Apply Burgess 2016 correction or MRlap -- MR-RAPS alone does not fix overlap-driven confounding (see Operational rule above), only weak-instrument bias |
| `cause()` runs forever | Default model fit on too many SNPs | Filter to sig SNPs (P < 1e-3) before `cause()`; `est_cause_params` uses the random subset |
| MAF column missing -> harmonise action=2 silently downgrades | EAF unavailable | Provide EAF or use `action = 3` and document the loss |

## Tool Installation Notes

`install.packages()` and `remotes::install_github()` lines for TwoSampleMR, MendelianRandomization, MVMR, MRPRESSO, mr.raps (CRAN-archived), cause, mrclust, lhcMR, DRMR, MRlap, simex, plus the plink2/1KG local-clumping note: `references/tool-installation.md`.

## STROBE-MR Reporting

STROBE-MR (Skrivankova 2021 JAMA 326:1614; BMJ 375:n2233): 20-item checklist required by Eur J Epi / Nat Genet / JAMA / Diabetologia since 2022. See causal-genomics/pleiotropy-detection for the per-item table -- that skill is the consolidated reporting reference for the full MR + sensitivity battery.

## References

Citations for every method named above: `references/bibliography.md`.

## Related Skills

- causal-genomics/colocalization-analysis - Confirm shared causal variant for cis-MR drug-target work
- causal-genomics/pleiotropy-detection - Deep dive on MR-PRESSO, Egger, contamination mixture diagnostics
- causal-genomics/fine-mapping - Credible-set construction at instrument loci before cis-MR
- causal-genomics/mediation-analysis - Two-step MR and MVMR difference method for X -> M -> Y mediation
- causal-genomics/transcriptome-wide-association - TWAS as MR-adjacent framework for gene-level inference
- causal-genomics/proteome-mr-drug-target - Dedicated cis-pQTL drug-target MR workflow with UKB-PPP/deCODE/Fenland and coloc triangulation
- population-genetics/association-testing - GWAS source for instrument selection
- clinical-databases/clinvar-lookup - Annotate instrument SNPs for downstream interpretation
