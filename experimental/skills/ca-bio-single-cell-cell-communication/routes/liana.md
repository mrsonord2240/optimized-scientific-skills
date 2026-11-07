# LIANA consensus
~~~sh
python scripts/liana_rank.py annotated.h5ad liana_results.csv --groupby cell_type --n-perms 1000
~~~
Provide an annotated object with log-normalized expression in X and at least two cell types. The source uses consensus resource, expr_prop=0.1 and use_raw=False; it reports magnitude and specificity ranks separately. Check current APIs and gene/resource compatibility; compare multiple resources and inspect abundance, depth and ambient RNA. Scores are co-expression proxies.
**Done when:** the full result, source/resource version, settings and biological limitations are retained.
