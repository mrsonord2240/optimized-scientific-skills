# CellPhoneDB specificity
~~~sh
python scripts/cellphonedb_statistical.py cellphonedb.zip meta.tsv counts_normalized.h5ad --out-dir cpdb_out --iterations 1000 --threads 1 --seed 1337
~~~
Provide a matching human resource, barcode metadata and normalized, unscaled counts. Source defaults use a 0.1 expression threshold and 0.05 p-value; permutation tests specificity, not direct binding. Source observed single-thread seed reproducibility only; verify installed behavior. See [CellPhoneDB notes](../references/cellphonedb.md).
**Done when:** output tables and resource, input, seed and limitations are recorded.
