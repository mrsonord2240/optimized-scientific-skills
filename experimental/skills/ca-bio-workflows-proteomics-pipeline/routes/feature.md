# Feature-level and DIA routes
~~~bash
Rscript scripts/msstats_maxquant.R evidence.txt proteinGroups.txt annotation.csv msstats_comparison.csv
Rscript scripts/diann_matrix.R report.parquet diann_log2_matrix.csv
~~~
Choose the matching route and a new output directory; scripts may overwrite. Done when the appropriate feature-level comparison or filtered log2 matrix is inspected. MSstats script hard-codes Control and Treatment; adapt the contrast matrix for other designs. DIA-NN produces a matrix, not a differential test; PG.MaxLFQ is already normalized. When MBR is enabled add the source-described Lib q-value filters: the script does not implement them. See [MSstats](../references/msstats.md) and [DIA-NN](../references/dia-nn.md).
