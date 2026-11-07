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
2. **Standard battery.** IVW (random-effects when Cochran Q p < 0.05) + MR-Egger (with NOME I^2_GX check) + weighted median + weighted mode + MR-PRESSO (NbDistribution `>= 10000` for stringent reporting). Report all five with point estimate, SE, p, 95% CI, and n_SNPs_used. Compute Egger I^2_GX; apply SIMEX if I^2_GX < 0.9 (see ../routes/simex.md).
3. **CHP escalation.** Trigger when rg > 0.3 OR PRESSO global p < 0.05 with > 50% nominal outliers OR Egger / median / mode disagree by > 2 SE. Run CAUSE (if `>= 100` significant SNPs after pruning) or LHC-MR (any N; uses genome-wide sumstats). Report ELPD delta + z + q (CHP fraction) + gamma (CHP-adjusted causal estimate).
4. **Triangulate.** Pre-MR Steiger filter; bidirectional MR (../routes/other.md); LCV gcp; LDSC rg report. Consensus across methods provides convergent evidence, subject to independent validation. Disagreement requires narrowing the scope (e.g., subgroup, cis-MR, time-varying analysis) rather than reporting a single point estimate.

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
| LDSC rg(exposure, outcome) `>= 0.3` or strong shared-factor biology | CAUSE (if sig SNPs `>= 100`; cause.md) OR LHC-MR (lhc-mr.md) | LCV gcp (lcv.md); cross-check IVW after Steiger filter |
| Many weak IVs (mean F < 20) | MR-RAPS with overdispersion + robust Huber loss (mr-raps-mr-clust.md) | MR-Mix or contamination mixture; report F-stat range |
| Suspected heterogeneous causal mechanisms (e.g. LDL on CHD via multiple lipoprotein pathways) | MR-Clust; report per-cluster IVW (mr-raps-mr-clust.md) | Pathway annotation of cluster instruments; Bayesian mixture |
| Cis-MR drug target (single locus, few SNPs in LD) | Colocalization (causal-genomics/colocalization-analysis) + Steiger | PWCoCo; conditional analysis; not Egger (low SNP count) |
| Polygenic exposure (heritability spread genome-wide; few significant loci) | LHC-MR (uses all SNPs; lhc-mr.md) | LDSC rg; genome-wide IVW with weak-IV-aware methods (RAPS) |
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

LCV gcp interpretation thresholds (0, 0.5, 0.6, 1) are tabulated in lcv.md.
