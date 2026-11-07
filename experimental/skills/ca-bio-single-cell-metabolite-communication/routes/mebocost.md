# MEBOCOST enzyme-sensor inference
~~~sh
python scripts/mebocost_run.py annotated.h5ad mebocost.conf results.csv --groupby cell_type --species human
~~~
Input X must already be log-normalized with gene symbols and annotated groups; verify normalization rather than inferring it from the maximum value. The config must reference locally available species-matched database files. Source settings use 50 cells per group for preprocessing, cutoff_prop 0.15, 1000 shuffles, seed 12345 and permutation FDR 0.05. The script prints excluded groups and writes the returned result. Check API and database version; retain resource provenance.
**Done when:** output, group sizes, filtering, species and permutation settings have been reviewed as hypotheses.
