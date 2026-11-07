# QC and impute
~~~bash
Rscript scripts/qc_cli.R extracted.rds qc.rds
~~~
Done when QC RSD falls while biological RSD is unchanged, wiped samples and missingness are reported, and the imputed matrix is NA-free. Source uses detection fraction 0.5, QCRSC minQC 5, QC RSD 30, PQN and seeded log-scale QRILC. It drops all-NA samples. D-ratio and blank checks are not implemented in this script and must be done separately; check [checkpoints](../references/checkpoints.md). The output sample_class is resynchronized; other metadata must be reconciled after dropped samples.
