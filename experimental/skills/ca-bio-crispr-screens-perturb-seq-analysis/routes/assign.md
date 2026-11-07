# Assign guides
~~~bash
python scripts/assign_sgrna.py sgrnas.h5ad --layer sgrna_counts --threshold 10 --out assigned.h5ad
~~~
Input layer is cells by guides with guide names as var_names. The source threshold is 10 reads; cells with zero or multiple qualifying guides are labelled none or multiplet. Inspect assignment rate, doublets and sparse/dense representation.
Done when assignments and ambiguous cells are reviewed.
