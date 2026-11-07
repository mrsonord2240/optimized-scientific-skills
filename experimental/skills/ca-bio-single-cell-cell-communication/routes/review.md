# Conditions, resource sensitivity and alternatives
The source condition-stability script permutes condition labels across cells within cell types. It is an exploratory cell-level sampling diagnostic, not a valid donor-level test when cells share donors or condition is assigned by donor. Use independent biological units and an appropriate replicate-aware comparison for condition inference; inspect scores and resource sensitivity rather than interpreting threshold crossings alone.
~~~sh
python scripts/condition_stability.py annotated.h5ad --cond-a control --cond-b treated --groupby cell_type --condition condition --n-null 10 --out stability.tsv
~~~
NicheNet requires a reliable receiver gene set and model assets; spatial methods, Tensor-cell2cell and Scriabin need separate workflows. See [NicheNet](../references/nichenet.md) and [methods](../references/methods.md). Do not use corrected expression as a substitute for raw-count, donor-aware differential analysis.
**Done when:** independent replication, resource sensitivity, confounds and orthogonal validation are documented.
Related shelf skills: bio-single-cell-cell-annotation, bio-single-cell-clustering and bio-single-cell-preprocessing.
