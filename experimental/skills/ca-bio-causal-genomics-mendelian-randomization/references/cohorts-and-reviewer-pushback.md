# Cohort gotchas and reviewer pushback

Read when picking a biobank/cohort for MR or answering a reviewer's MR critique. Verbatim from SKILL.md "Cohort Gotchas" and "Anticipated Reviewer Pushback".

## Cohort Gotchas

- UKB GWAS commonly include non-EUR participants; pan-UKB EUR/AFR/EAS/SAS subsets are separate releases. Mixing ancestries inflates instrument strength via population stratification rather than biology.
- FinnGen DF12 (2024) cohort: Finnish founder effects produce narrower LD blocks and higher winner's curse magnitude than UKB at matched sample size.
- MVP / GBMI / AoU: multi-ancestry meta-analyses; stratify by ancestry before MR or use ancestry-specific subsets.
- UKB-on-UKB MR creates one-sample-equivalent bias regardless of "different GWAS file" appearance; use MRlap or move the outcome to an external cohort (FinnGen, BBJ, MVP).

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Sample overlap between exposure and outcome GWAS?" | LDSC bivariate intercept reported; MRlap or sample-overlap-corrected IVW applied; document overlap fraction |
| "Weak instruments (F < 10)?" | F computed from EXPOSURE; per-instrument and mean F reported; MR-RAPS used as sensitivity if mean F borderline |
| "Horizontal pleiotropy?" | IVW + Egger + weighted median + weighted mode + MR-PRESSO; if rg > 0.3 also CAUSE (see pleiotropy-detection) |
| "Reverse causation?" | Steiger filter applied; bidirectional MR ran; LCV gcp reported if rg > 0.3 |
| "Pre-registered?" | OSF protocol filed; STROBE-MR all 20 items reported |
| "InSIDE assumption?" | INstrument Strength Independent of Direct Effect -- pleiotropic effects alpha uncorrelated with instrument-exposure effects gamma. Tested via Egger intercept + CHP-aware sensitivity (CAUSE) |
