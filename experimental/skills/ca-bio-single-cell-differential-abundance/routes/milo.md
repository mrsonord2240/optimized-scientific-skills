# Milo neighborhoods
~~~sh
Rscript scripts/milo_run.R input_sce.rds results.csv
~~~
Input SingleCellExperiment needs a suitable PCA reduced dimension and sample, condition and cell_type metadata. Each sample must map to one condition and at least two independent samples per condition are required by this wrapper; more may be needed for stable inference. Source settings: seed 42, k=30, d=30, prop=0.1 and refined neighborhoods. Inspect graph, sample counts and SpatialFDR. Add study-specific covariates only after checking model rank and design.
**Done when:** neighborhood results, sample-level design, settings and limitations are recorded.
