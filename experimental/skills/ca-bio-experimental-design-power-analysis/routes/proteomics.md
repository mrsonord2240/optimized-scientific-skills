# Proteomics multiplicity approximation
~~~sh
Rscript scripts/proteomics_power.R scenarios.csv power.csv
~~~
CSV columns are d, power, alpha and n_proteins. The script calculates per-group sample size for a two-sample t-test using raw alpha and the source's Bonferroni alpha/number-of-proteins comparison. Gaussian model, effect size, independence assumptions where relevant, missingness and the family definition require review; this is not an FDR simulation.
**Done when:** both rounded sample-size estimates and assumptions are retained.
