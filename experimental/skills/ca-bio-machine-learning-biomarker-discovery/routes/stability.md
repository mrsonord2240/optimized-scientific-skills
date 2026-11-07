# Subsampled selection stability
~~~sh
python scripts/stability_run.py samples.csv outcome frequencies.csv stability.json
~~~
CSV rows must be independent samples with numeric features and binary outcome. The source helper scripts/stability_selection.py standardizes inside each half-sample and uses a fixed feature budget derived from pi=0.6 and expected-error setting 1; its theoretical bound has assumptions that may fail in correlated omics. The wrapper reports real and one permuted-label control; repeat permutations and validate on independent data. Stability alone does not prove signal or causality.
**Done when:** frequencies, stability index, null comparison and assumptions are reviewed.
