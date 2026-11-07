# Run low-MOI SCEPTRE
~~~bash
python scripts/run_sceptre_safe.py input.rds results.tsv
~~~
Input RDS supplies response and guide matrices, target metadata, covariates and discovery pairs; see [SCEPTRE notes](../references/sceptre-low-moi.md). Use a fresh output path; the wrapper refuses an existing result and requires a successful worker plus a valid newly produced table. Inspect results.tsv.run.json and worker_returncode; successful completion and structural validation do not establish statistical calibration. Source worker scripts/run_sceptre.R is invoked by the wrapper.
Done when a fresh output, worker status, covariates, calibration and missing values are reviewed.
