# Load 10x Cell Ranger output

Python:

```bash
python scripts/load_10x.py raw_feature_bc_matrix.h5 loaded.h5ad
python scripts/load_10x.py filtered_feature_bc_matrix/ loaded.h5ad var_names=gene_ids
```

R:

```bash
Rscript scripts/load_10x.R filtered_feature_bc_matrix/ loaded.rds
```

`load_10x.py` takes a matrix directory or an `.h5`, keeps non-GEX features (`gex_only=False`, antibody capture and CRISPR guides land in `var['feature_types']`), makes names unique, prints cells x genes and the feature types, and writes a gzip-compressed `.h5ad`. `var_names=gene_ids` (default) is Ensembl IDs, stable across annotation releases; `gene_symbols` is ambiguous and release-dependent. For `.h5` input the names are symbols and the IDs are in `var['gene_ids']`; `var_names` is ignored. `load_10x.R` returns a Seurat object with `min.cells=3`, `min.features=200`; with several feature types it keeps `Gene Expression` as RNA and adds the others as extra assays.

- Raw versus filtered: store `raw_feature_bc_matrix` for ambient-RNA tools; `filtered_feature_bc_matrix` is cell-called only.
- `make_unique` appends `-1`/`-2` to duplicate symbols and can mask distinct paralog or PAR loci; ID names avoid that.
- `genome='GRCh38'` selects one genome from a multi-species `.h5` in `sc.read_10x_h5`.
- `sc.read_10x_mtx(cache=True)` caches by path only: rereading with a different `var_names` returns the stale object. Delete the cache or omit it.
- Multi-sample work: load each sample separately and keep QC per sample (preprocessing Skill) before merging.

Done when cells x genes and the feature types are reported and the output file exists.
