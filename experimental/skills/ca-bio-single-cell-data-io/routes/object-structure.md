# Object structure

## AnnData

`X` and every `layers[...]` are cells x genes. `obs` and `obsm` annotate cells; `var` and `varm` annotate genes; `obsp` and `varp` are square pairwise graphs; `uns` is unstructured.

```python
adata.layers['counts'] = adata.X.copy()   # integer UMIs, for scVI, DESeq2 and HVG seurat_v3
# ... normalize_total + log1p fill X ...
adata.raw = adata                         # frozen log-normalized full-gene snapshot, survives later var-subsetting
adata = adata[:, adata.var['highly_variable']].copy()
```

- `adata.raw` has its own `var`; `adata.raw.to_adata()` rebuilds it. Slicing by obs also slices raw; slicing by var does not shrink raw.
- Tools disagree on `X` versus `.raw.X` (`use_raw=`). Set it explicitly and confirm which holds counts and which holds log values.
- Check sparsity with `scipy.sparse.issparse(adata.X)`.
- `read_h5ad(backed='r')`: only `X` edits persist; `obs`, `var` and `obsm` edits are not written. Write to a new file.

## Seurat v5

Counts, data and scale.data are layers of an `Assay5`; v4 used fixed slots (`GetAssayData(slot=)`).

```r
counts <- LayerData(seurat_obj, layer = 'counts')   # v5; GetAssayData(slot=) is the superseded v4 form
counts <- seurat_obj[['RNA']]$counts                # shorthand
merged <- merge(obj1, y = c(obj2, obj3), add.cell.ids = c('S1', 'S2', 'S3'))
merged <- JoinLayers(merged)                        # merge() splits layers (counts.1, counts.2); rejoin first
```

Rejoin after `merge()` before any function that expects one layer.

Done when each matrix is in its named slot and the layer or `use_raw` choice is stated.
