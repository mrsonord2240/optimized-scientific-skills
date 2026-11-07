# Sensitivity battery
~~~sh
Rscript scripts/sensitivity_battery.R harmonized.tsv battery.rds
~~~
Supply a TwoSampleMR harmonized table including the beta, standard-error, allele, sample-size and frequency fields required by the tests. Retained instruments are selected by mr_keep when present. The workflow computes IVW, MR-Egger, weighted median and weighted mode estimates, heterogeneity, Egger intercept, leave-one-out and Steiger diagnostics. Steiger requires adequate sample size and allele-frequency data. Interpret these as complementary sensitivity analyses, not proof that pleiotropy is absent. Review the [method decision framework](../references/methods-and-thresholds.md). See [failure modes](../references/failure-modes.md) and [reporting](../references/reporting.md).
**Done when:** the RDS contains estimates and diagnostics, and instrument selection and assumptions have been reviewed.
