# Bidirectional MR and Steiger directionality

Read when testing causal direction (reverse MR, `steiger_filtering()`, `directionality_test()`) or when a Steiger flag contradicts prior biology. Verbatim from SKILL.md "Bidirectional and Steiger" and "Steiger filter false flag under unmeasured confounding".

## Bidirectional and Steiger

```r
exposure_rev <- format_data(outcome_raw, type = 'exposure')  # treat former outcome as exposure
outcome_rev <- format_data(exposure_raw, type = 'outcome')
dat_rev <- harmonise_data(exposure_rev, outcome_rev, action = 2)
results_rev <- mr(dat_rev, method_list = 'mr_ivw')

dat_filt <- steiger_filtering(dat)  # per-SNP; flags SNPs where variance(Y) > variance(X)
dir_test <- directionality_test(dat) # global; correct_causal_direction == TRUE if forward
if (is.null(dir_test)) stop("directionality_test() returned NULL -- dat needs samplesize_col set ",
                             "on both read_*_data()/format_data() calls upstream; see the Standard ",
                             "Workflow section above.")
```

**Report direction operationally:** forward p < 5e-8 with reverse p > 0.05; `directionality_test()` `correct_causal_direction == TRUE` with Steiger p < 0.05; point estimate in reverse direction has |effect| substantially smaller than forward (reflecting reverse-instrument-strength asymmetry). Example sentence: "Forward MR showed BMI -> T2D (IVW beta = 0.85, p = 2e-15); reverse MR was null (IVW beta = 0.02, p = 0.61); Steiger directionality test favored the forward direction (p = 3e-7)."

### Steiger filter false flag under unmeasured confounding

**Trigger:** Applying `steiger_filtering()` on traits with unmeasured shared confounders (e.g. SES).

**Mechanism:** Steiger compares variance explained in exposure vs outcome per SNP; an unmeasured confounder upstream of both produces SNPs that explain more variance in the outcome than the exposure, falsely flagging "reverse causation" (Lutz 2022 Genet Epidemiol 46:139).

**Symptom:** Many SNPs flagged as wrong-direction yet biology and prior MR support forward causation.

**Fix:** Treat Steiger as a heuristic, not gospel; cross-validate direction with bidirectional MR (forward + reverse with independent instrument sets); for known-confounder-rich domains (psychiatric traits, SES proxies) use LCV or LHC-MR instead, which jointly model confounders.
