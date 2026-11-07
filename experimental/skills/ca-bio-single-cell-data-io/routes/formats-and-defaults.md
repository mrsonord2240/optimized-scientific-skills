# Storage formats and surprising defaults

| Format | Use when | Weak when |
|--------|----------|-----------|
| h5ad (HDF5) | Default single-machine Python I/O and sharing | Not cloud-native; partial and concurrent reads limited |
| zarr | Cloud or object store, larger than memory, lazy (Dask), anndata 0.11+ v3 sharding | Many small files on a local filesystem; v2/v3 version skew breaks old readers |
| RDS | Seurat-only workflow, full object fidelity | R only; version-tied |
| h5mu (MuData) | Multimodal (RNA + ADT + ATAC), one AnnData per modality | Fewer tools read it; call `mdata.update()` |
| Loom | Legacy interchange (velocyto, older Seurat) | `write_loom(write_obsm_varm=False)` drops obsm/varm by default; aging |

Defaults that change the data:

| Call | Surprise | Do |
|------|----------|----|
| `sc.read_10x_mtx(gex_only=True)` | Drops Antibody, CRISPR and custom features | `gex_only=False`, then split on `var['feature_types']` |
| `sc.read_10x_mtx(var_names='gene_symbols')` | Non-unique, release-dependent | `'gene_ids'` for joins |
| `AnnData.write_h5ad(compression=None)` | No compression | `compression='gzip'` |
| `write_loom(write_obsm_varm=False)` | Embeddings lost | Pass `True` |
| `read_h5ad(backed='r')` | Only `X` edits persist | Write to a new file |
| sceasy `drop_single_values=TRUE` | Constant columns deleted | Set `FALSE` |
| `sc.read_10x_mtx(cache=True)` | Cache keyed by path only | Delete the cache when `var_names` changes |

Multimodal objects (CITE-seq, Multiome): see the multimodal-integration Skill. Spatial data: spatial-data-io.

Done when the format choice and any default overridden are stated.
