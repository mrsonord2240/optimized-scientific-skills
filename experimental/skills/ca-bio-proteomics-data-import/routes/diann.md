# Import DIA-NN protein matrix
~~~bash
python scripts/load_diann.py report.parquet log2_matrix.tsv
~~~
The source uses Q.Value, PG.Q.Value and Global.PG.Q.Value <=0.01 then pivots PG.MaxLFQ by protein group and run with first-value aggregation, converts zero to NaN and log2-transforms. Check report schema and duplicate group/run values; first is not a reconciliation method. For library-based MBR the source notes a Lib.PG.Q.Value alternative. Version and downstream missingness require review.
Done when q-filters, input schema and matrix are reviewed.
