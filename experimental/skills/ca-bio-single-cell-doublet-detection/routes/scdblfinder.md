# scDblFinder
~~~sh
Rscript scripts/scdblfinder_run.R input_sce.rds capture_id scored.rds 20260924
~~~
Input is a SingleCellExperiment with raw counts and a capture-level metadata column. The script uses per-capture scoring and a seeded serial BiocParallel parameter; source defaults infer rate from cell count. If samples were multiplexed in one capture, preserve capture identity and account for total lane loading rather than treating each demultiplexed subset as a physical capture. Check installed API and loading assumptions.
**Done when:** score, class, capture counts, rate assumption and seed are retained and suspect populations reviewed.
