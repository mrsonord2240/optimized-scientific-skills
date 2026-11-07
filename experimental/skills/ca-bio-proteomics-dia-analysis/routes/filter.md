# Filter an existing DIA-NN report
~~~bash
python scripts/filter_report.py diann_out/report.parquet log2_matrix.tsv --context global
~~~
The source filters run-level precursor and protein-group plus global precursor and protein-group q-values at <=0.01 before pivoting PG.MaxLFQ by Protein.Group and Run. For the source's 1.9.x MBR case use --context lib, requiring Lib.Q.Value and Lib.PG.Q.Value. Zero becomes NaN before log2; duplicate group/run rows use first aggregation, so inspect disagreements. Report-versus-matrix row counts do not establish matrix FDR.
Done when schema, q-level/context, duplicate behavior and matrix are reviewed.
