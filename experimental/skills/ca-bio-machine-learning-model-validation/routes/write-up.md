# What to state in the write-up

No script. State each item; mark an unavailable one as unavailable with a reason.

1. Cohort, outcome and horizon, prevalence, intended use.
2. Split design: unit of independence, folds, repeats, seeds, any external or temporal test.
3. Every preprocessing, selection, tuning, calibration and threshold step, and where each was fit.
4. Discrimination with an uncertainty interval; for tuned models the nested estimate beside the flat score and the number of configurations searched.
5. Calibration intercept, slope, reliability curve and Brier on untouched data.
6. Net benefit at pre-specified thresholds when a decision is claimed.
7. Subgroup performance with subgroup sizes and uncertainty.
8. Missing data, exclusions, limitations, code and data versions.

Reporting standard for a clinical prediction claim: TRIPOD+AI.

- Cross-validation on one source estimates performance on new patients from that source. A different time, place or setting is external validation; calibration degrades first. With several cohorts, use leave-one-site-out (`routes/site-or-time.md`).
- Apparent performance overstates the future; a development calibration slope below 1 is the optimism signal. Shrink the model or penalise.
- Sample size: do not use "10 events per variable". Size with the Riley minimum-sample-size framework (`pmsampsize`). When p is much larger than n, those formulas do not apply: use heavy penalisation and nested validation.
- Choosing the features themselves, or a treatment-effect estimate, is outside this Skill.

Done when items 1 to 8 are each stated or marked unavailable.
