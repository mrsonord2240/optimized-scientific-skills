# Scanpy marker ranking
~~~sh
python scripts/scanpy_markers.py clustered.h5ad markers.csv --groupby leiden
~~~
Input X must contain suitable log-normalized expression and at least two groups. The source uses Wilcoxon with expressing fractions; results are descriptive and depend on the cluster-vs-rest composition. Check effect size and in/out fractions, gene identifiers and installed API before interpreting. Do not use cell-level cluster tests as donor-level condition DE.
**Done when:** full rankings, input layer, test and marker specificity have been reviewed.
