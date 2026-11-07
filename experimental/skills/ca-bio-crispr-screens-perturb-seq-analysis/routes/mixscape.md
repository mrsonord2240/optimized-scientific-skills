# Classify escapers
~~~bash
python scripts/mixscape_filter.py normalized.h5ad --pert-key sgrna_assignment --control NTC --out ko_cells.h5ad
~~~
Use normalized log1p RNA and the actual control label; default neighbors 20 and random_state 0. The script writes KO cells only, excluding both NP and controls. Preserve the original control cells separately and combine appropriately for downstream comparisons; do not use KO-only output as if it includes controls. Weak effects may be misclassified.
Done when control identity, KO retention and weak-phenotype limitations are reviewed.
