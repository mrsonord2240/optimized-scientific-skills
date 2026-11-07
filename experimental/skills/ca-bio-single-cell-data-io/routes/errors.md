# Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Converted object has genes and cells swapped | Transpose missing, or applied without swapping the metadata axis | Transpose the matrix and move `obs` to column metadata and `var` to row metadata |
| Layers, embeddings or `raw` missing after conversion | Lossy converter dropped non-`X` slots | Diff inventories (`scripts/h5ad_inventory.py`); `.raw` via `schard::h5ad2sce(use.raw=TRUE)`; reattach the rest |
| Cannot run EmptyDrops, SoupX or CellBender | Only the filtered matrix was kept | Re-obtain the RAW Cell Ranger matrix |
| ADT or guide counts absent after loading | `gex_only=True` default | Reload with `gex_only=False`, split by `var['feature_types']` |
| Session dies reading a large object | Dense materialization of a sparse matrix | Stay sparse; `backed='r'` (Python) or on-disk layers (Seurat v5) |
| Downstream tool uses the wrong values | Read `X` where `.raw.X` was expected | Set `use_raw=` explicitly |
| Duplicate gene symbols suffixed `-1`/`-2` | `make_unique` on distinct loci | Load with `var_names='gene_ids'` |

Done when the symptom is gone and the slot inventory matches the original.
