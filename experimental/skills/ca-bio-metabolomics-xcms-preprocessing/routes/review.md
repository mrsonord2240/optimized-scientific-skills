# Review and downstream QC

```bash
Rscript examples/xcms_workflow.R
```

This source demonstration uses bundled faahKO data and deletes its temporary table; it is not a user-data route. For actual files use [workflow](workflow.md) and inspect EICs and alignment. The source describes matched-filter detection, QC-anchored alignment, CAMERA grouping and XCMS QC filters, but these require input-specific adaptation and are not included in the file wrapper. Treat filled values as model-derived, preserve unfilled values and compute per-feature fill fraction. CAMERA grouping can both merge and split compounds; do not equate feature counts with compound counts. Use `bio-metabolomics-normalization-qc` for drift and QC workflows.

Done when parameter choices, filled fractions, redundancy and QC limitations are reviewed before downstream interpretation.
