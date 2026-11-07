# Mixscape response signatures
~~~sh
python scripts/mixscape_run.py annotated.h5ad mixscape.h5ad
~~~
Input needs normalized expression in X and perturbation and gene_target metadata, with NT controls. The source computes a local signature using 20 neighbors and batch_size 1000, then fits Mixscape classes. This script densifies sparse X as in the source; verify sufficient memory before running. An NP call can reflect ineffective perturbation or no detectable response; do not equate it with gene function. Inspect fractions per target and preserve raw counts separately for donor-aware downstream inference.
**Done when:** the object, class fractions, controls and limitations are reviewed.
