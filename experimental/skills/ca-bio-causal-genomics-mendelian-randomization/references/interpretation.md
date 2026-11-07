# Reconciliation: When Methods Disagree
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
