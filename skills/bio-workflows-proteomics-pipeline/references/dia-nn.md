# DIA-NN workflow

Read this when the search output is DIA-NN `report.parquet`: q-value filtering before pivoting, the PG.MaxLFQ matrix, and the 0 -> NA step before log2.

### DIA-NN Workflow
DIA-NN 1.9+ defaults to report.parquet (the only default in 2.0); read it with arrow, not read.delim. Filter on q-values BEFORE pivoting, or low-confidence rows enter the matrix. Route to proteomics/dia-analysis for the mechanics.
Run `scripts/diann_matrix.R report.parquet diann_log2_matrix.csv` (q-value filters, `PG.MaxLFQ` pivot, 0 -> NA, log2). `PG.MaxLFQ` is already normalized, so go straight to limma with no re-normalization.
