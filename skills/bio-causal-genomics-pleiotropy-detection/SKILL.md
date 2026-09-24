---
name: bio-causal-genomics-pleiotropy-detection
description: Detect and adjust for horizontal pleiotropy in two-sample Mendelian randomization by distinguishing uncorrelated (UHP) from correlated (CHP) pleiotropy and choosing among Egger, MR-PRESSO, MR-RAPS, CAUSE, LHC-MR, LCV, MR-Clust, MR-Mix, and contamination-mixture methods. Use when validating an MR causal claim, running the STROBE-MR sensitivity battery, suspecting a shared heritable confounder, working under weak-instrument or polygenic-exposure regimes, or reconciling discordant estimates across robust methods.
tool_type: r
primary_tool: TwoSampleMR
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: TwoSampleMR 0.5.11+, MendelianRandomization 0.9.0+, MR-PRESSO 1.0+, CAUSE 1.2.0+, MR-Clust 0.1.0+, MRMix 0.1+, mr.raps 0.4.1+ (GitHub), LHC-MR 0.0.0.9000+ (GitHub), LCV (script-based, no version tag), simex 1.8+.

```r
install.packages(c('remotes', 'TwoSampleMR', 'MendelianRandomization', 'simex'))
remotes::install_github('rondolab/MR-PRESSO')          # CRAN-never; GitHub-only
remotes::install_github('jean997/cause')               # CHP-aware (Morrison 2020)
remotes::install_github('cnfoley/mrclust')             # Mechanism heterogeneity
remotes::install_github('gqi/MRMix')                   # Mixture-of-distributions
remotes::install_github('qingyuanzhao/mr.raps')        # CRAN-archived 2025-03; install from GitHub
remotes::install_github('LizaDarrous/lhcMR')           # Heritable-confounder MR
# LCV: git clone https://github.com/lukejoconnor/LCV (R scripts, no package)
```

Inputs are typically a harmonized TwoSampleMR data.frame (beta.exposure, beta.outcome, se.exposure, se.outcome, SNP, effect_allele, eaf, etc.). LHC-MR and LCV take genome-wide GWAS sumstats with LDSC-style merged SNPs; CAUSE takes pruned signature SNPs plus full sumstats for nuisance estimation.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- For GitHub-only packages, check the repo HEAD vs the local install date

If code throws errors, introspect the installed package and adapt the example rather than retrying. For reproducibility, pin TwoSampleMR, MRPRESSO, and CAUSE versions in the methods section of any report, and record harmonization choices.

# Pleiotropy Detection in Mendelian Randomization

**"Validate my MR result against pleiotropic bias"** -> Decompose violations of the exclusion-restriction assumption into uncorrelated horizontal pleiotropy (UHP, addressable by Egger / median / mode / MR-PRESSO) and correlated horizontal pleiotropy (CHP, addressable only by CAUSE / LHC-MR / LCV), then run a method battery whose assumptions span both regimes.

- R: `TwoSampleMR::mr()` (IVW + Egger + median + mode), `mr_pleiotropy_test()`, `mr_heterogeneity()`, `mr_leaveoneout()`, `directionality_test()`
- R: `MRPRESSO::mr_presso()` for UHP outlier removal + distortion test
- R: `cause::cause()` for CHP-aware estimation; `mrclust::mr_clust_em()` for mechanism-heterogeneous instruments
- R: `MendelianRandomization::mr_conmix()` for contamination mixture; `MRMix::MRMix()` for mixture-of-distributions

**Practice boundary:** Every estimate this Skill produces is a population-average causal effect from GWAS summary statistics. It is research/epidemiological output, not a clinical tool. If a request turns a population MR result into an individual treatment or diagnostic decision ("should I personally start drug X"), decline the personal recommendation, name the population-vs-individual gap explicitly (baseline risk, comorbidities, and effect-modifier interactions a GWAS-derived average cannot capture), redirect to a physician, and offer to continue the population-level analysis instead.

## UHP vs CHP: The Central Postdoc-Grade Distinction

Horizontal pleiotropy comes in two regimes, and most "standard" MR sensitivity methods address only one of them.

| Regime | Definition | InSIDE assumption | Methods that handle it |
|--------|------------|--------------------|------------------------|
| UHP (uncorrelated horizontal pleiotropy) | Pleiotropic effect alpha_j independent of instrument-exposure effect gamma_j | Holds | IVW (balanced UHP only), MR-Egger, weighted median, weighted mode, MR-PRESSO, MR-RAPS, MR-Mix, contamination mixture |
| CHP (correlated horizontal pleiotropy) | alpha_j correlates with gamma_j through a shared upstream factor (heritable confounder, network mediator) | Violated | CAUSE, LHC-MR, LCV, MR-Clust (partial), Steiger-filtered MR (partial) |

**InSIDE = INstrument Strength Independent of Direct Effect** (Bowden 2015 IJE 44:512). Plain English: across SNPs, the per-SNP pleiotropic effect alpha and per-SNP instrument-exposure effect gamma are treated as independent random variables. CHP is the case where they covary because both flow from a shared upstream genetic factor.

**The trap (Morrison 2020 Nat Genet 52:740):** IVW, MR-Egger, MR-PRESSO, and GSMR are all blind to CHP. Under a shared heritable confounder they each return a plausible-looking corrected causal estimate that is systematically biased in the direction of the confounder. The MR-PRESSO global test does not flag CHP because correlated pleiotropy is not an outlier pattern, it is a population mean shift in the alpha distribution conditional on gamma.

**Operational rule:** If genetic correlation rg(exposure, outcome) is high (LDSC `>= 0.3`) or biology strongly suggests a shared upstream factor, the IVW / Egger / PRESSO triple is insufficient. Add CAUSE (preferred when sig SNPs `>= 100`) or LHC-MR (preferred for polygenic genome-wide IVs).

## Operational Decision Flow (4 Steps)

1. **Compute genetic correlation (LDSC).** Run `ldsc.py --rg <exposure.sumstats.gz>,<outcome.sumstats.gz>` (see causal-genomics/genetic-correlation). If `|rg| > 0.3`, CHP is plausible -> flag for Step 3 escalation. If the LDSC rg standard error spans zero broadly, treat low-rg evidence as weak rather than confirming absence of CHP.
2. **Standard battery.** IVW (random-effects when Cochran Q p < 0.05) + MR-Egger (with NOME I^2_GX check) + weighted median + weighted mode + MR-PRESSO (NbDistribution `>= 10000` for stringent reporting). Report all five with point estimate, SE, p, 95% CI, and n_SNPs_used. Compute Egger I^2_GX; apply SIMEX if I^2_GX < 0.9 (see examples/simex_egger_correction.R).
3. **CHP escalation.** Trigger when rg > 0.3 OR PRESSO global p < 0.05 with > 50% nominal outliers OR Egger / median / mode disagree by > 2 SE. Run CAUSE (if `>= 100` significant SNPs after pruning) or LHC-MR (any N; uses genome-wide sumstats). Report ELPD delta + z + q (CHP fraction) + gamma (CHP-adjusted causal estimate).
4. **Triangulate.** Pre-MR Steiger filter; bidirectional MR (examples/bidirectional_mr.R); LCV gcp; LDSC rg report. Consensus across methods supports a publication-ready claim. Disagreement requires narrowing the scope (e.g., subgroup, cis-MR, time-varying analysis) rather than reporting a single point estimate.

## Algorithmic Taxonomy

| Method | Models | UHP-robust | CHP-robust | Min #SNPs | Fails when | Citation |
|--------|--------|-----------|-----------|-----------|------------|----------|
| Inverse-variance weighted (IVW) | Weighted regression through origin | Balanced UHP only | No | 3 | Directional UHP; CHP; weak IV bias; heterogeneity | Burgess 2013 Genet Epidemiol 37:658 |
| MR-Egger intercept + slope | IVW + free intercept | Directional UHP | No | >=10 for power | NOME violated (I^2_GX < 0.9); <10 SNPs; CHP | Bowden 2015 IJE 44:512 |
| Weighted median | Median of Wald ratios | Up to 50% invalid | No | >=10 | >50% invalid; CHP | Bowden 2016 Genet Epidemiol 40:304 |
| Weighted mode (MBE) | Mode of estimate density | Plurality valid | Partial | >=10 | Multimodal estimates from CHP clusters | Hartwig 2017 IJE 46:1985 |
| Cochran Q | Heterogeneity across Wald ratios | Total heterogeneity flag, not direction-specific | No | 3 | Cannot distinguish UHP from heterogeneity from CHP | Del Greco M F 2015 Stat Med 34:2926 |
| MR-PRESSO | Detect + remove UHP outliers via RSS-out | Yes (assumes majority valid) | No | >=4 | >50% pleiotropic; any CHP; small n | Verbanck 2018 Nat Genet 50:693 |
| GSMR + HEIDI-outlier | Outlier removal via single-instrument estimate heterogeneity | Yes | No | >=10 | CHP (HEIDI-outlier is heterogeneity-driven) | Zhu 2018 Nat Commun 9:224 |
| MR-RAPS | Profile likelihood with overdispersion + Huber/Tukey loss | Yes; weak-IV robust down to F~10 | Partial via overdispersion | >=10 | Strong CHP; at extreme-weak IV (mean F well below 10, e.g. ~2-3) the overdispersion estimator can degenerate ("very small"/negative, forced to tau2=0) and RAPS can underperform plain IVW -- do not assume RAPS beats IVW without checking its warnings | Zhao 2020 Ann Stat 48:1742 |
| MR-Mix | Mixture-of-distributions over valid + invalid | Yes | Partial | >=20 | Few SNPs; very heterogeneous CHP | Qi & Chatterjee 2019 Nat Commun 10:1941 |
| Contamination mixture | Profile likelihood over contamination fraction | Yes | Partial | >=20 | Few SNPs | Burgess 2020 Nat Commun 11:376 |
| MR-Clust | k-means over Wald estimates with NULL cluster | Yes | Diagnostic for CHP via clusters | >=20 | Single-mechanism exposure (no clustering signal) | Foley 2021 Bioinformatics 37:531 |
| CAUSE | Bayesian mixture: shared causal + shared-factor (CHP) components | Yes | Yes (explicit) | >=100 sig SNPs at p<5e-8 | <100 sig SNPs; non-overlapping GWAS samples | Morrison 2020 Nat Genet 52:740 |
| LHC-MR | Latent heritable confounder + bidirectional + heritability | Yes | Yes | Genome-wide GWAS sumstats | Heritability mis-estimated; severe sample overlap | Darrous 2021 Nat Commun 12:7274 |
| LCV (latent causal variable) | gcp parameter on genome-wide rg | N/A (not an MR method) | Diagnostic | Genome-wide GWAS sumstats | Heritability low; non-Gaussian effect distribution | O'Connor & Price 2018 Nat Genet 50:1728 |

Methodology evolves; verify against the Burgess & Thompson textbook (2nd ed 2021), Hemani 2018 (basic four-method battery), and Sanderson 2022 Nat Rev Methods Primers 2:6 before locking a sensitivity battery.

## Decision Tree by Scenario

| Scenario | Primary estimator | Sensitivity / triangulation |
|----------|-------------------|----------------------------|
| Many strong IVs, no biological shared trait suspected | IVW + Egger + weighted median + weighted mode + PRESSO | Cochran Q; leave-one-out; F-stat; Steiger filtering |
| LDSC rg(exposure, outcome) `>= 0.3` or strong shared-factor biology | CAUSE (if sig SNPs `>= 100`; references/cause.md) OR LHC-MR (references/lhc-mr.md) | LCV gcp (references/lcv.md); cross-check IVW after Steiger filter |
| Many weak IVs (mean F < 20) | MR-RAPS with overdispersion + robust Huber loss (references/mr-raps-mr-clust.md) | MR-Mix or contamination mixture; report F-stat range |
| Suspected heterogeneous causal mechanisms (e.g. LDL on CHD via multiple lipoprotein pathways) | MR-Clust; report per-cluster IVW (references/mr-raps-mr-clust.md) | Pathway annotation of cluster instruments; Bayesian mixture |
| Cis-MR drug target (single locus, few SNPs in LD) | Colocalization (causal-genomics/colocalization-analysis) + Steiger | PWCoCo; conditional analysis; not Egger (low SNP count) |
| Polygenic exposure (heritability spread genome-wide; few significant loci) | LHC-MR (uses all SNPs; references/lhc-mr.md) | LDSC rg; genome-wide IVW with weak-IV-aware methods (RAPS) |
| Reverse causation suspected | Bidirectional MR with Steiger filter; LHC-MR (jointly estimates both directions) | directionality_test; effect-size r2 comparison |
| Population-level summary discordant with biology | Re-examine instrument selection; check Winner's curse; LD pruning settings | Triangulate with cis-MR; family-based MR if available |

## Quantitative Thresholds

| Metric | Threshold | Source / rationale |
|--------|-----------|--------------------|
| F-statistic per SNP | `>=10` strong; `<10` weak | Burgess 2011 IJE 40:755 (rule of thumb); weak-IV bias toward observational confounded estimate |
| I^2_GX (NOME) | `>=0.9` Egger reliable | Bowden 2016 IJE 45:1961 |
| I^2_GX (NOME) intermediate | 0.6-0.9 SIMEX-corrected Egger | Bowden 2016 IJE |
| I^2_GX (NOME) severe | `<0.6` drop Egger; use MR-RAPS or CAUSE | Bowden 2016 IJE; SIMEX unreliable below 0.6 |
| Egger min SNP count for adequate power | `>=10` | Bowden 2015 IJE 44:512 |
| Cochran Q significance | p < 0.05 indicates heterogeneity | Del Greco M F 2015 Stat Med 34:2926 |
| MR-PRESSO NbDistribution | 1000 exploratory; `>=5000` publication; `>=10000` stringent | Verbanck 2018 Nat Genet 50:693 (Methods) |
| MR-PRESSO global test p | < 0.05 -> heterogeneity / outliers present | Verbanck 2018 Nat Genet |
| MR-PRESSO distortion test p | < 0.05 -> outliers materially shifted estimate; if `>= 0.05` report uncorrected IVW | Verbanck 2018 Nat Genet |
| MR-PRESSO min instruments | `>=4` to run; `>=10` for non-degenerate global test | Verbanck 2018 Nat Genet |
| MR-PRESSO SignifThreshold | 0.05 default | Verbanck 2018 Nat Genet |
| MR-PRESSO majority-valid breakdown | Fails when `>50%` instruments pleiotropic | Verbanck 2018 Nat Genet (theoretical limit) |
| Weighted median validity | Robust to `<=50%` invalid IVs | Bowden 2016 Genet Epidemiol 40:304 |
| Weighted mode validity | Plurality-valid (largest valid subset is most common estimate) | Hartwig 2017 IJE 46:1985 |
| CAUSE min #SNPs | `>=100` p < 5e-8 SNPs after pruning | Morrison 2020 Nat Genet 52:740 (Supplement) |
| CAUSE delta_ELPD criterion | one-sided p < 0.05; z = delta_elpd / se(delta_elpd); z > 1.96 standard; z > 3.0 stringent | Morrison 2020 Nat Genet |
| LDSC rg suggesting CHP | `>= 0.3` flags need for CAUSE / LHC-MR | Operational rule; see causal-genomics/genetic-correlation |
| Steiger r^2 difference | Reverse-causal flag at any per-SNP r2_GY > r2_GX | Hemani 2017 PLoS Genet 13:e1007081 |
| Standard sensitivity battery | IVW + Egger + median + mode + PRESSO + Steiger + LOO | Hemani 2018 eLife 7:e34408 / STROBE-MR 2021 |

LCV gcp interpretation thresholds (0, 0.5, 0.6, 1) are tabulated in references/lcv.md.

## Standard Sensitivity Battery (Working Reference)

**Goal:** Run the canonical UHP-focused MR sensitivity suite on harmonized two-sample data.

**Approach:** Compute IVW + Egger + median + mode side-by-side; test Egger intercept and heterogeneity; run MR-PRESSO with `>=5000` distributions for publication or `>=10000` for stringent reporting; apply Steiger filter; leave-one-out; report all estimates.

Run it from `examples/sensitivity_battery.R` (IVW / Egger / weighted median / weighted mode via `mr(dat, method_list = ...)`, `mr_heterogeneity`, `mr_pleiotropy_test`, `Isq`, `directionality_test`, `mr_presso`, `mr_leaveoneout`, contamination mixture, MR-RAPS, STROBE-MR table); it simulates its own `dat`, so replace that block with your harmonized data.frame. Keep `set.seed(42)` ahead of the whole battery: `mr()`'s weighted-median and weighted-mode bootstrap SEs and `mr_presso()`'s global/outlier tests are all Monte-Carlo, so an unseeded rerun changes those SEs and p-values. The example uses `NbDistribution = 5000`; raise it to 10000 for stringent reporting. SIMEX correction is separate: examples/simex_egger_correction.R.

## Bidirectional MR Procedure

1. **Forward MR:** instrument exposure E, test effect on outcome Y (primary)
2. **Reverse MR:** instrument outcome Y, test effect on exposure E (using outcome-direction instruments)
3. **Steiger pre-filter both directions:** `steiger_filtering(dat)`; drop SNPs where outcome r^2 > exposure r^2 before primary IVW
4. **Compare estimates:** null reverse + significant forward strengthens the forward causal claim; bidirectional significance flags feedback / shared confounder / reciprocal causation
5. **LCV gcp orthogonal check:** genome-wide directional inference independent of the instrument set

Working code: examples/bidirectional_mr.R.

**Interpretation cheat-sheet:**

| Forward p | Reverse p | Reading |
|-----------|-----------|---------|
| significant | non-significant | Forward causal claim strengthened |
| non-significant | significant | Re-examine instrument-exposure assignment; the "outcome" may causally drive the "exposure" |
| significant | significant | Feedback loop, shared confounder, or reciprocal causation; resolve with LHC-MR |
| non-significant | non-significant | No evidence of causation in either direction |

When forward and reverse both clear Steiger and both IVW p < 0.05, run LHC-MR jointly rather than reporting two univariable estimates.

## Reconciliation Across Methods

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| IVW and Egger agree (small Egger intercept); median and mode agree | Likely true causal; minimal pleiotropy | Report all; STROBE-MR; emphasize agreement |
| IVW significant; Egger non-significant with similar slope; PRESSO global p > 0.05 | Egger underpowered (few SNPs) OR Egger NOME violated | Check I^2_GX; SIMEX-correct if 0.6 < I^2_GX < 0.9 |
| IVW shifted relative to Egger / median / mode; Egger intercept significant | Directional UHP | Trust Egger slope; PRESSO-corrected IVW; mode estimator |
| IVW + Egger + median + mode + PRESSO all agree but CAUSE delta_ELPD non-significant or in opposite direction | CHP via shared confounder | Trust CAUSE; report all five UHP methods alongside but flag the discordance; check LDSC rg |
| All UHP methods agree; LCV gcp ~ 0 | Genetic correlation only, not causation | Causation evidence is weak; report rg explicitly; consider colocalization (cis-MR) instead |
| MR-Clust shows >=2 distinct non-null clusters | Heterogeneous mechanisms | Report per-cluster estimates; do not summarize as a single effect |
| Steiger fails on a substantial fraction of instruments | Reverse causation OR exposure measurement error | Run bidirectional MR; check exposure GWAS heritability; LHC-MR |
| MR-PRESSO global p < 0.05 but corrected estimate similar to uncorrected | >50% pleiotropic OR CHP masquerading as UHP | Re-prune instruments; switch to weighted-mode / CAUSE / LHC-MR |

**Operational rule for publication:** Report IVW (primary), Egger slope + intercept, weighted median, weighted mode, MR-PRESSO global + distortion + corrected, Cochran Q, Steiger directionality, F-statistic distribution, I^2_GX (for Egger validity), and at least one CHP-aware method (CAUSE or LHC-MR) when rg `>= 0.3` or biology suggests shared upstream. Failure to report a CHP-aware result when CHP is plausible is a reviewer-flagged red flag since 2020.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| MR-PRESSO crashes with `Not enough intrumental variables` | Fewer than 4 SNPs | Need >=4 for PRESSO; for cis-MR with few SNPs use colocalization |
| Egger intercept p < 0.05 but I^2_GX = 0.5 | NOME violated; intercept is artifactually inflated | SIMEX-correct or do not trust Egger; use MR-RAPS instead |
| `Isq()` not found | TwoSampleMR version where Isq is unexported | Compute manually: Q_GX = sum((beta_GX/se_GX)^2); I2 = (Q_GX - (n-1))/Q_GX, clipped to [0,1] |
| MR-RAPS `package not found` after CRAN install | CRAN-archived 2025-03-01 | `remotes::install_github('qingyuanzhao/mr.raps')`; call via `TwoSampleMR::mr_raps()` wrapper or `mr.raps::mr.raps()` directly (MendelianRandomization does NOT export `mr_raps`) |
| MR-RAPS warns `overdispersion parameter is very small` / `negative, using tau2 = 0` and underperforms IVW | Extreme-weak IV (mean F well below the 10 threshold); the profile-likelihood overdispersion fit degenerates | Report both RAPS and IVW with the warning text; do not switch to RAPS as primary at extreme-weak F without checking which one is closer to other robust estimates (median/mode) |
| CAUSE delta_ELPD CI spans zero; Pareto-k > 0.7 | <100 sig SNPs OR severe sample overlap | Use LHC-MR; or report CAUSE with the explicit caveat |
| CAUSE crashes `Error in -1 * comp[2, 1] : non-numeric argument to binary operator` | `loo` >= 2.6 changed `loo_compare()`'s output layout; `cause`'s internal `in_sample_elpd_loo()` still reads the pre-2.6 layout by position (see references/cause.md) | Pin `loo` to <2.6 (e.g. 2.5.1) in a library resolved ahead of `cause`'s |
| Steiger labels most instruments reverse-causal | Exposure GWAS imprecise OR sample size mismatch | Treat as one signal; cross-check with bidirectional MR |
| LHC-MR runtime > 24h | Default nCores=1 on full sumstats | Use nCores >= 4; restrict to LDSC-overlapping SNPs first |
| MR-PRESSO outliers all on same chromosome | Genome-wide LD not properly pruned; clumping window too narrow | Re-clump at r^2 < 0.001 in 10 Mb window |
| MR-Mix / contamination mixture return a confident point estimate below their documented 20-SNP minimum (not NA) | Small n destabilizes the mixture-model fit without an internal check that refuses to report | Check n before trusting the output: below 20 SNPs, treat MR-Mix / conmix point estimates and CIs as unreliable regardless of apparent significance; report n_SNPs alongside the estimate; cross-check against MR-Clust and the full UHP battery rather than the mixture methods alone |

## Reference Files

| File | Read when |
|------|-----------|
| references/failure-modes.md | A method's result looks wrong or unexpected: MR-Egger NOME / SIMEX (incl. the SIMEX overshoot caveat), PRESSO majority-outlier breakdown or CHP false negative, Egger with few SNPs, Steiger inversion, CAUSE with few SNPs, LCV under non-Gaussian effects |
| references/cause.md | Running CAUSE, or interpreting q / eta / gamma / delta_ELPD / Pareto-k |
| references/mr-raps-mr-clust.md | Weak instruments (MR-RAPS loss function, overdispersion, nested `parameters =` calling form) or heterogeneous mechanisms (MR-Clust) |
| references/lhc-mr.md | Polygenic exposure, severe sample overlap or bidirectional interest: LHC-MR workflow and the CAUSE-vs-LHC-MR choice table |
| references/lcv.md | LCV gcp: `RunLCV` usage (field names) and the gcp interpretation table |
| references/reporting.md | Writing up: instrument and sensitivity-battery supplementary tables, anticipated reviewer pushback, STROBE-MR items |

## References

- Bowden J et al 2015 Int J Epidemiol 44:512 (MR-Egger; InSIDE)
- Bowden J et al 2016 Int J Epidemiol 45:1961 (NOME, I^2_GX, SIMEX)
- Cook JR & Stefanski LA 1994 JASA 89:1314 (SIMEX framework)
- Burgess S 2020 Nat Commun 11:376 (contamination mixture)
- Burgess S & Thompson SG 2021 (Mendelian Randomization 2nd ed)
- Darrous L et al 2021 Nat Commun 12:7274 (LHC-MR)
- Foley CN et al 2021 Bioinformatics 37:531 (MR-Clust)
- Hartwig FP et al 2017 Int J Epidemiol 46:1985 (weighted mode)
- Hemani G et al 2017 PLoS Genet 13:e1007081 (Steiger orientation)
- Hemani G et al 2018 eLife 7:e34408 (TwoSampleMR framework)
- Morrison J et al 2020 Nat Genet 52:740 (CAUSE; CHP)
- O'Connor LJ & Price AL 2018 Nat Genet 50:1728 (LCV)
- Sanderson E et al 2022 Nat Rev Methods Primers 2:6 (MR Primer)
- Skrivankova VW et al 2021 JAMA 326:1614 (STROBE-MR statement); BMJ 375:n2233 (explanation)
- Verbanck M et al 2018 Nat Genet 50:693 (MR-PRESSO)
- Zhao Q et al 2020 Ann Stat 48:1742 (MR-RAPS)

## Related Skills

- causal-genomics/mendelian-randomization - Primary causal estimation that this sensitivity battery validates
- causal-genomics/genetic-correlation - LDSC rg required for Step 1 of the decision flow; CHP escalation trigger
- causal-genomics/colocalization-analysis - Required for cis-MR drug-target signals where instruments are too few for Egger
- causal-genomics/fine-mapping - Identify causal variants underlying instrument loci
- causal-genomics/mediation-analysis - Multivariable MR for mediator-adjusted causal estimates
- population-genetics/association-testing - GWAS summary statistics underlying MR instruments
- clinical-biostatistics/effect-measures - Translate MR estimates to clinical effect measures
