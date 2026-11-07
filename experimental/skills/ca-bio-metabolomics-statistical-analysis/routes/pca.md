# PCA and Hotelling screening

```bash
Rscript scripts/pca_hotelling.R peaks.csv pareto 0.05 t2.csv
```

CSV has feature IDs in the first column, features in rows and samples in columns, with no missing values. The source fits ropls PCA with explicit Pareto scaling and computes an F-based Hotelling T2 threshold; inspect component count and assumptions, and do not automatically delete flagged samples. Compare scaling choices and inspect batch and injection order. Output may be overwritten. Tight pooled-QC clustering alone does not validate the whole run.

Done when preprocessing, scaling, components, flags and independent QC evidence are documented.
