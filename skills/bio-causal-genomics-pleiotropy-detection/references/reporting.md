# Reporting: Supplementary Tables, Reviewer Pushback, STROBE-MR

## Required Supplementary Tables

**Instrument table** (one row per SNP retained for primary analysis):

| Column | Content |
|--------|---------|
| rsID, chr, pos | Variant identifier and genome position |
| EA, OA, EAF | Effect allele, other allele, effect allele frequency in exposure GWAS |
| beta_E, se_E, p_E, F | Exposure-side estimate, SE, p-value, per-SNP F-statistic |
| beta_Y, se_Y, p_Y | Outcome-side estimate, SE, p-value (harmonized to EA) |
| harmonization_action | 1=kept, 2=flipped, 3=dropped (palindromic ambiguity) |
| palindromic_flag | TRUE / FALSE; tracked per Hartwig 2016 IJE 45:1717 |
| Steiger_direction | forward / reverse / inconclusive |
| Steiger_p | per-SNP directionality p-value |

**Sensitivity-battery table** (one row per method):

| Column | Content |
|--------|---------|
| method | IVW / Egger / WM / WMode / PRESSO (raw + corrected) / RAPS / CAUSE |
| estimate, se, p, 95% CI | Point estimate and inference |
| n_SNPs_used | Post-harmonization, post-Steiger SNP count |
| heterogeneity_p | Q for IVW; Q' for Egger; global p for PRESSO |
| intercept_p | Egger only (directional UHP test) |
| ELPD_delta + z | CAUSE only (sharing - causal; negative + |z| > 1.96 -> causal) |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Was CHP checked for?" | LDSC rg reported (causal-genomics/genetic-correlation); if rg > 0.3, CAUSE or LHC-MR ran; q posterior reported |
| "Why CAUSE and not LHC-MR?" | CAUSE preferred when `>= 100` significant SNPs available (Morrison 2020). LHC-MR preferred when significant-SNP set is small or polygenic, using genome-wide sumstats (Darrous 2021) |
| "Egger NOME?" | I^2_GX computed; if 0.6 <= I^2_GX < 0.9, SIMEX correction applied; if < 0.6, Egger dropped in favor of MR-RAPS |
| "PRESSO doesn't catch CHP?" | Confirmed (Morrison 2020); CAUSE / LHC-MR reported alongside PRESSO for that reason |
| "Steiger filter applied pre-MR or post?" | Pre-MR: SNPs failing per-SNP Steiger directionality dropped before primary IVW |
| "Why no replication cohort?" | Two-sample design uses independent exposure and outcome cohorts; if same biobank, MRlap used or noted as a limitation |
| "Why not just trust the IVW?" | IVW assumes balanced UHP and no CHP; both violated routinely; sensitivity battery is the standard since STROBE-MR 2021 |
| "Effect size is implausibly large" | Re-examine F-statistic distribution for weak IV bias; check Winner's curse; consider Wald ratio at a single strong instrument as sanity check |

## STROBE-MR Reporting (Skrivankova 2021)

| Item | Required content |
|------|-----------------|
| 1-3 | Title / abstract / background indicates this is an MR study; pre-registered protocol |
| 4-7 | Study design, data sources, instrument selection criteria (p-threshold, LD clumping, MAF) |
| 8-11 | Harmonization, palindrome handling, allele alignment |
| 12-14 | F-statistic distribution; weak-instrument bias mitigation |
| 15-17 | Primary MR method + all sensitivity methods + CHP-aware method when relevant |
| 18-19 | Pleiotropy tests, Steiger, heterogeneity |
| 20 | Discussion of remaining assumption violations; limitations |

Sub-items (30 total) detail per-method reporting. The full statement (JAMA 326:1614) and explanation (BMJ 375:n2233) are now reviewer-required at most cardiovascular and psychiatric journals since 2022.
