# Assess representation and replicates
~~~bash
python scripts/library_representation.py screen.count.txt --out library_representation.tsv
~~~
Use a MAGeCK TSV with guide ID first, Gene and numeric sample columns. Validate nonnegative numeric values, unique guide IDs, no missing data and identify all-zero samples before using this helper; it does not enforce all source checks. It reports zero and low counts, quantiles and p90/p10 skew; zero p10 yields NaN. For Gini, within-condition correlations and plots run scripts/screen_qc.py after editing its fixed stage_map and condition_map and placing screen.count.txt in the working directory; adjust script path as needed. Stage-specific limits are in [metrics](../references/metrics.md).
Done when sample mappings, exclusions, stage thresholds and outliers are reviewed.
