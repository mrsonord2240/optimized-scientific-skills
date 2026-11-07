# SingleR with a local reference
~~~sh
Rscript scripts/singler_local.R query.rds reference.rds label.main result.rds
~~~
Query and reference must be compatible SingleCellExperiment objects with appropriate normalized expression and matching feature identifiers; the reference needs the named label column. The script saves labels, pruned labels and scores. Choose an appropriate reference and review pruning and score distributions; confident nearest-label assignment does not establish that the true type is present.
**Done when:** predictions and uncertainty are saved and checked against reference scope, markers and QC.
