# Estimate surrogate variables
~~~sh
Rscript scripts/sva_from_files.R expression.csv metadata.csv sv.rds
~~~
Expression CSV has feature IDs in its first column and samples in remaining columns; metadata CSV has id and condition. Sample names and metadata IDs must match exactly. The helper removes features with any non-finite value, reports removal, builds condition and intercept models, and uses the source's leek estimator. Complete-feature filtering can bias toward abundant features; alternatively use justified imputation before analysis. Include surrogate variables in the inference model rather than treating a cleaned matrix as raw observations.
**Done when:** the RDS, feature exclusion count, model and downstream interpretation are reviewed.
