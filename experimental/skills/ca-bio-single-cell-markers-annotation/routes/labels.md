# Apply reviewed labels
~~~sh
python scripts/apply_labels.py clustered.h5ad labels.csv annotated.h5ad --groupby leiden
~~~
CSV requires cluster and label columns. Labels must be curated for this tissue, reference and dataset; unmatched clusters become Unassigned. Avoid assigning cell types from fixed illustrative cluster numbers. Confirm with multiple marker patterns and flag uncertainty.
**Done when:** labels, unmapped clusters and supporting evidence are saved.
