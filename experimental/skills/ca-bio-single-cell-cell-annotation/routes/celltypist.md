# CellTypist
~~~sh
python scripts/celltypist_annotation.py --input clustered.h5ad --model model.pkl --output annotated.h5ad --figure annotation.png --counts-output counts.csv
~~~
Input must contain raw counts in layers/counts, seeded leiden labels and X_umap; the script normalizes counts to CP10K and log1p and uses majority voting. The model must already exist locally, with species, tissue and gene identifiers matched. The source example marks confidence above 0.5 as high; that is illustrative, not a universal calibrated threshold. Inspect distributions and verify installed APIs.
**Done when:** labels, confidence, reference version, counts and ambiguous populations are reviewed against markers and QC.
