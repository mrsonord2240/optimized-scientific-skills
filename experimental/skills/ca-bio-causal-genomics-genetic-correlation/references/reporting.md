# Reconciliation Across Methods
| Pattern | Likely cause | Action |
|---------|--------------|--------|
| LDSC rg and HDL rg agree (independent samples) | Both methods converging on true value | Report HDL as primary (lower SE); LDSC as sensitivity |
| LDSC rg substantially below HDL rg, LDSC intercept large | Sample overlap; HDL is biased toward phenotypic correlation | Report LDSC as primary; flag overlap; do not report HDL |
| Global LDSC rg ~ 0 but LAVA shows multiple Bonferroni-significant local rg | Locus-level cancellation in global average | Report both; the biology is "shared at specific loci, divergent overall" |
| LAVA significant but univariate h2 non-significant at hit locus | Spurious bivariate without identified local h2 signal | Filter univariate first; do NOT report bivariate at unidentified loci |
| Popcorn rho_ge << 1 across many trait pairs | Population-specific causal architecture | Real finding; report rho_ge alongside within-pop h2 |
| LDSC ratio > 0.2 | Population stratification or model misfit | Re-check ancestry; consider LD-score reference mismatch; report with caveat |
| LCV gcp ~ 0 with large rg | Genetic correlation without (partial) causation in either direction | Shared confounder hypothesis is preferred; do NOT report as causal |
| LCV gcp > 0 (significant) with large rg | Partial-to-full causation in tested direction | Combine with bidirectional MR + CHP-aware sensitivity; this is supportive but not sufficient |

**Operational rule for publication:** Report LDSC rg + intercept as primary global statistic; report HDL only if overlap is verified < 5%; complement with LAVA local rg whenever global rg is near zero or biology suggests heterogeneity; report LCV gcp when downstream MR is planned; trans-ancestry analyses require Popcorn (not within-population LDSC).

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| LDSC `category not found` after `--rg` | LD score column header mismatch (custom reference) | Inspect M_5_50 and `.l2.ldscore` headers; align with `--ref-ld-chr` prefix |
| LDSC ratio > 1 (negative h2 z-score) | Severe stratification or wrong LD-score ancestry | Switch to ancestry-matched reference; check for population structure |
| HDL convergence warning / NA SE | Reference panel mismatch or extreme overlap | Verify SVD eigen reference path; switch to LDSC when overlap suspected |
| LAVA `Insufficient SNPs at locus` for most loci | LD reference and partition file from different builds | Match GRCh37 vs GRCh38; align LD reference to partition file |
| LAVA bivariate rg = +/- 1 at boundary | Univariate filter not applied; locus is unidentified | Apply `univ$p < 0.05/N_loci` filter to BOTH traits before `run.bivar()` |
| Popcorn complains about MAF format | sumstats EAF column missing or NA | Provide EAF; do not impute from external reference (creates miscalibration) |
| HESS `--num-shared` defaulting to wrong value | Forgot to set explicitly; default 0 is independent | Always set explicitly; if unknown, use LDSC intercept to infer overlap |
| GenomicSEM `ldsc()` returns negative-definite covariance | Numerical instability with many traits | Inspect per-pair LDSC results; drop low-h2 traits; regularize |
| rg point estimate > 1 with CI overlapping 1 | Sampling variance; same-trait pair near identity | Report as "rg not distinguishable from 1"; constrained likelihood at the rg=1 boundary gives different SE -- LRT against H0: rg=1 is more precise than Wald CI |
| Binary-vs-continuous trait pair scale concern | Reviewer asks about liability-vs-observed scale propagation | LDSC rg is scale-invariant -- case-control h2 liability vs observed scale propagates equivalently into rg; no correction needed |

## Required Reporting for rg Analyses

| Component | Required |
|-----------|----------|
| Per-trait h2 + SE + intercept + mean chi-square | Yes |
| Bivariate rg + SE + p | Yes |
| gcov_int (cross-trait intercept) | Yes; non-zero under known overlap is expected, not bias |
| LD reference panel + ancestry | Yes |
| Method used (LDSC / HDL / LAVA / Popcorn) | Yes; rationale per Decision Tree |
| Local rg supplementary (LAVA) | If global rg null but biology suggests sharing |
| Sample-size: Neff per trait | Yes |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Sample overlap?" | LDSC: gcov_int reported; non-zero under known overlap is expected, NOT bias. HDL: only valid if independent (<5% overlap) |
| "Why LDSC not HDL?" | HDL gives lower variance but is biased > 5% overlap; LDSC is the conservative default |
| "Local vs global rg?" | If global rg modest but biology suggests sharing, LAVA (Werme 2022) reported as supplementary |
| "Cross-ancestry?" | Popcorn for trans-ancestry; rg < 1 in trans is real biology, not noise |
| "Does rg motivate CHP-MR?" | If `|rg| > 0.3`, CAUSE / LHC-MR sensitivity reported (cross-ref pleiotropy-detection) |
| "rg = 1 boundary?" | If CI includes 1, reported as "not distinguishable from rg=1"; constrained LRT alternative provided |
