---
name: bio-single-cell-data-io
description: Read, write, create, and convert single-cell objects across AnnData (Python), Seurat (R), and SingleCellExperiment (R). Use when loading 10X Cell Ranger output (raw vs filtered), importing or exporting h5ad/RDS/h5mu/zarr, building AnnData or Seurat objects from matrices, moving objects between Python and R, or debugging lost layers, transposed matrices, or mangled gene names during conversion.
tool_type: mixed
primary_tool: Seurat
license: MIT
---

## Version Compatibility

Reference examples tested with: scanpy 1.10+, Seurat 5.0+, anndata 0.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Single-Cell Data I/O

**"Load my 10X data"** -> Parse a Cell Ranger matrix into an annotated object (cells, genes, counts, metadata).
- Python: `sc.read_10x_mtx()` / `sc.read_10x_h5()` -> AnnData
- R: `Read10X()` / `Read10X_h5()` -> `CreateSeuratObject()`

## Installation

```bash
pip install scanpy anndata
# multimodal: pip install muon mudata
```

```r
install.packages('Seurat')
# Conversion (prefer maintained tools; SeuratDisk is abandoned):
remotes::install_github('scverse/anndataR')          # pure-R h5ad/zarr I/O + conversion; requires R >= 4.5
BiocManager::install('zellkonverter')                # SCE <-> AnnData
remotes::install_github('cellgeni/schard')           # robust pure-R h5ad reading
```

## Governing Principle

The dominant failure in single-cell I/O is not a crash; it is a silent semantic change to the matrix during read or conversion. Three traps drive almost every lost-data bug.

Orientation is opposite across ecosystems. AnnData is cells x genes (observations on rows, `obs` indexes rows, `var` indexes columns); Seurat and SingleCellExperiment are genes x cells (features on rows, cells on columns). So `adata.X` is the transpose of `LayerData(seu)` and `assay(sce)`. A faithful conversion must transpose AND swap which axis the metadata annotates; getting the transpose right but the metadata axis wrong is the single most common silent conversion bug.

Sparse storage compounds the transpose. R `Matrix::dgCMatrix` is CSC; scanpy conventionally stores `X` as CSR. Transposing a CSR matrix yields CSC for free, so a correct AnnData->Seurat hop involves both a logical transpose and a CSR<->CSC change. Forcing dense (`.toarray()`, `as_dense=` on write) on a 500k-cell x 30k-gene float32 matrix materializes ~60 GB; keep `X` and layers sparse and check with `scipy.sparse.issparse(adata.X)`.

Conversion is lossy by default, and the loss is silent. Cross-ecosystem hops drop `layers`, `obsp`/`varp`, nested `uns`, and coerce categoricals to character/NA. `adata.raw` has its own `var` (its purpose is to survive HVG subsetting) and tools disagree on whether to read `X` or `.raw.X` (`use_raw=`), so a mismatched expectation silently uses the wrong matrix. Always diff slot inventories before and after any cross-ecosystem conversion, and keep the original file.

One more governing fact: the Cell Ranger filtered matrix is cell-CALLED, not ambient-corrected. The widespread claim that the filtered matrix is "decontaminated" is false. Keep the RAW (unfiltered) matrix, because EmptyDrops, SoupX, CellBender, and DecontX all require it and filtered-only storage is irreversible.

## Choosing a Storage Format

| Format | Backing | Use when | Fails / weak when |
|--------|---------|----------|-------------------|
| h5ad (HDF5) | single file | Default single-machine Python I/O and sharing | Not cloud-native; concurrent/partial reads limited |
| zarr (directory of chunks) | object store | Cloud/S3, larger-than-memory, parallel/lazy (Dask), anndata 0.11+ v3 sharding | Many small files awkward on local FS; v2/v3 version skew breaks old readers |
| RDS | single R binary | Seurat-only workflow, full object fidelity in R | R-only; not portable to Python; version-tied |
| h5mu (MuData) | HDF5 | Multimodal (RNA + ADT + ATAC), one AnnData per modality | Less tool support than h5ad; needs `mdata.update()` discipline |
| Loom (HDF5) | single file | Legacy interchange (velocyto, older Seurat) | `write_loom(write_obsm_varm=False)` DROPS obsm/varm by default; aging |

Methods and tool maturity move fast here. Before committing a conversion route, verify the chosen package is still maintained and matches installed versions (`packageVersion`, `pip show`).

## Loading 10X Cell Ranger Output

**Goal:** Read a Cell Ranger matrix correctly, keeping the raw matrix and non-GEX features when present.

**Approach:** Read the raw (unfiltered) MEX/HDF5 matrix; select stable Ensembl IDs for reproducible joins; retain Antibody/CRISPR features by disabling `gex_only`.

```python
import scanpy as sc

# raw_feature_bc_matrix has every barcode (needed by EmptyDrops/SoupX/CellBender); filtered_feature_bc_matrix has only called cells
adata = sc.read_10x_mtx('raw_feature_bc_matrix/', var_names='gene_ids', gex_only=False)
# gene_ids (Ensembl) is stable across annotation releases; gene_symbols (default) is ambiguous and non-unique
# gex_only=False keeps Antibody Capture / CRISPR Guide; split later by adata.var['feature_types']
adata.var_names_make_unique()
```

```r
library(Seurat)
counts <- Read10X(data.dir = 'filtered_feature_bc_matrix/')          # list when multiple feature types present
seurat_obj <- CreateSeuratObject(counts = counts, project = 'PBMC', min.cells = 3, min.features = 200)
```

Read functions return symbols by default. `sc.read_10x_h5` has no `var_names` argument (symbols by default; Ensembl IDs land in `var['gene_ids']`). `make_unique` appends `-1`/`-2` to duplicate symbols, which can mask distinct paralog/PAR loci, so prefer IDs when joining datasets.

## AnnData Object Structure

**Goal:** Place counts, normalized values, metadata, and embeddings in the conventional slots so downstream tools find them.

**Approach:** Keep integer counts in `layers['counts']`, log-normalized values in `X`, and a frozen full-gene snapshot in `.raw` before HVG subsetting.

```python
import anndata as ad

# X is (n_obs, n_vars) = cells x genes; obs indexes rows, var indexes columns
adata.layers['counts'] = adata.X.copy()   # integer UMIs, kept to recompute or feed count models (scVI, DESeq2)
# ... normalize_total + log1p populate X ...
adata.raw = adata                         # frozen log-normalized full-gene snapshot; survives later var-subsetting
adata = adata[:, adata.var['highly_variable']].copy()
```

Slot roles: `X`/`layers` align to both axes (each exactly cells x genes); `obs`/`obsm` align to cells; `var`/`varm` align to genes; `obsp`/`varp` are square pairwise graphs; `uns` is unstructured. `adata.raw.to_adata()` reconstitutes the snapshot. Slicing the parent by obs also slices raw on obs, but var-slicing does NOT shrink raw.

## Seurat v5 Object Structure

**Goal:** Read and write the right assay layer under the v5 layers API.

**Approach:** Use `LayerData()`/`$`-accessors; rejoin split layers after `merge()` before any function expecting one layer.

```r
counts <- LayerData(seurat_obj, layer = 'counts')      # v5; GetAssayData(slot=) is the superseded v4 form
counts <- seurat_obj[['RNA']]$counts                   # shorthand
merged <- merge(obj1, y = c(obj2, obj3), add.cell.ids = c('S1', 'S2', 'S3'))
merged <- JoinLayers(merged)                           # merge() splits layers (counts.1, counts.2); rejoin first
```

Seurat v5 stores `counts`/`data`/`scale.data` as layers in an `Assay5`; v4 used fixed slots via `GetAssayData(slot=)`. After `merge()`, layers split per object until `JoinLayers()`.

## Converting Between Python and R

**Goal:** Move an object across ecosystems without dropping layers, embeddings, or `raw`.

**Approach:** Prefer a maintained pure-R or Python-pinned converter; transpose and remap metadata; diff slots before and after.

| Tool | Direction | Maintained 2026 | Use when |
|------|-----------|-----------------|----------|
| anndataR | AnnData <-> SCE <-> Seurat; h5ad+zarr R/W | Yes (v1.2.0, pure R, no Python) | First choice for R-native, Python-free h5ad/zarr I/O and conversion. Requires R >= 4.5 (`getRversion()`); falls back to zellkonverter/schard on R 4.4 |
| zellkonverter | AnnData <-> SCE | Yes (Bioc 3.23) | Mature SCE<->AnnData; robust Python reader with pinned anndata. `raw=TRUE` can silently fail to populate `altExp` on h5ad written by current anndata (0.13.x) -- verify with `length(altExpNames(sce)) > 0`; if empty, use schard's `use.raw=TRUE` below |
| schard | h5ad -> Seurat/SCE (read-only) | Yes | Robust pure-R READING of h5ad (SeuratDisk replacement); the reliable route to recover `.raw` when zellkonverter's `raw=TRUE` fails (`use.raw=TRUE`). Reads `X` only -- does not carry a `counts` layer; pair with zellkonverter or reattach the counts matrix manually if one is needed downstream |
| anndata2ri | AnnData <-> SCE (rpy2) | Yes | Live mixed Python+R sessions / Jupyter `%%R` |
| sceasy | everything -> AnnData hub | Aging | Quick one-call conversion (mind `drop_single_values` data loss) |
| SeuratDisk | AnnData <-> h5Seurat | NO (last commit 2023, broken on Seurat v5) | Avoid for new work; legacy only |

```r
# Preferred R-native read of an h5ad written in Python (no reticulate)
library(anndataR)
adata <- read_h5ad('data.h5ad')
seurat_obj <- adata$to_Seurat()
# Or via Bioconductor with a pinned Python anndata:
# library(zellkonverter); sce <- readH5AD('data.h5ad'); writeH5AD(sce, 'out.h5ad')
```

**Recovering `.raw`.** `zellkonverter::readH5AD(reader='R', raw=TRUE)` is the documented route but is unreliable: verified against a real h5ad with a genuine `.raw` group (anndata 0.13.x write, zellkonverter 1.16.0), it returns an SCE with `altExpNames` empty -- the `.raw` group exists in the file but the R-side reader does not surface it. `schard::h5ad2sce(path, use.raw=TRUE)` reliably recovers it instead:

```r
library(schard)
sce_raw <- schard::h5ad2sce('data.h5ad', use.raw = TRUE)   # full-gene raw snapshot
sce_hvg <- schard::h5ad2sce('data.h5ad', use.raw = FALSE)  # default: X as stored (e.g. HVG-subsetted)
length(SingleCellExperiment::altExpNames(sce_hvg))          # runnable "diff slot inventories" check for zellkonverter's raw=TRUE -- 0 means it silently failed
```

**Seurat -> h5ad.** No tool converts this direction directly; the route below is `as.SingleCellExperiment()` + `zellkonverter::writeH5AD()`, verified end to end on a real Seurat 5.5.0 PBMC object (h5ad reloaded correctly in Python):

```r
library(Seurat)
library(zellkonverter)

sce <- as.SingleCellExperiment(seurat_obj)   # counts -> counts, data -> logcounts, reductions -> reducedDims
writeH5AD(sce, 'out.h5ad')
```

`scale.data` is silently dropped by `as.SingleCellExperiment()`: a Seurat object with `counts`/`data`/`scale.data` layers keeps only `counts`/`logcounts` (`assayNames(sce)`), and the resulting h5ad has no trace of `scale.data` after reload in Python. Reductions (e.g. PCA) survive into `obsm`. Save `scale.data` separately first if it is needed downstream: `as.matrix(LayerData(seurat_obj, layer='scale.data'))`.

zellkonverter maps asymmetrically: `obsm`->`reducedDims`, `varm`->a `rowData` matrix column (NOT reducedDims), `obsp`/`varp`->`colPairs`/`rowPairs`, `uns`->`metadata()` (lossy). sceasy's `drop_single_values=TRUE` silently deletes every obs/var column with one unique value (a one-sample object loses its constant batch/condition label), so set `FALSE`.

## API Defaults That Surprise

| Call | Surprising default | Consequence |
|------|--------------------|-------------|
| `sc.read_10x_mtx(gex_only=True)` | drops Antibody/CRISPR/Custom features | CITE-seq ADT and guides silently vanish; set `gex_only=False` |
| `sc.read_10x_mtx(var_names='gene_symbols')` | non-unique, release-dependent symbols | Use `'gene_ids'` for reproducible cross-dataset joins |
| `AnnData.write_h5ad(compression=None)` | no compression (gzip default removed after v0.6.16) | Larger files; pass `compression='gzip'` |
| `write_loom(write_obsm_varm=False)` | obsm/varm dropped | Embeddings lost on Loom write |
| `read_h5ad(backed='r')` | only `X` edits persist | `obs`/`var`/`obsm` edits in backed mode are NOT written; re-`.write()` to a new file |
| `sceasy convertFormat(drop_single_values=TRUE)` | constant columns deleted | Single-value batch/condition labels lost; set `FALSE` |
| `sc.read_10x_mtx(cache=True)` | cache keyed by path only | Re-reading a path with different `var_names` returns the STALE object; delete the `.h5ad` cache or omit `cache` |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Converted object has genes and cells swapped | Transpose not applied (or applied without swapping metadata axis) | Transpose the matrix AND move `obs`<->col-meta, `var`<->row-meta |
| Layers / embeddings / `raw` missing after conversion | Lossy converter dropped non-`X` slots | Diff slot inventories; for `.raw` specifically use `schard::h5ad2sce(use.raw=TRUE)`, not zellkonverter's `raw=TRUE` (unreliable, see Converting section); re-attach other slots manually |
| Cannot run EmptyDrops/SoupX/CellBender | Only the filtered matrix was kept | Re-obtain and store the RAW (unfiltered) Cell Ranger matrix |
| ADT/guide counts absent after loading 10X | `gex_only=True` (default) dropped non-GEX features | Reload with `gex_only=False`, split by `var['feature_types']` |
| Kernel/session dies reading a large object | Dense materialization of a sparse matrix | Keep sparse; use `backed='r'` (Python) or BPCells/on-disk layers (Seurat v5) |
| Downstream tool uses wrong values | Tool read `X` vs `.raw.X` against expectation | Set `use_raw=` explicitly; confirm which matrix holds counts vs lognorm |
| Duplicate gene symbols collapsed or suffixed oddly | `make_unique` appended `-1`/`-2` to distinct loci | Load with `var_names='gene_ids'` for stable identifiers |

## Related Skills

- single-cell/preprocessing - QC, normalization, and HVG selection after loading
- single-cell/doublet-detection - per-sample doublet calling on raw counts after loading
- single-cell/clustering - dimensionality reduction and clustering on the loaded object
- single-cell/multimodal-integration - MuData/h5mu handling for CITE-seq and Multiome
- spatial-transcriptomics/spatial-data-io - SpatialData/zarr I/O for spatial omics
- workflows/scrnaseq-pipeline - end-to-end scRNA-seq pipeline that starts from data loading

## References

- Virshup I, et al. (2023) The scverse project provides a computational ecosystem for single-cell omics. Nature Biotechnology 41:604-606. DOI 10.1038/s41587-023-01733-8
- Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA (2024) anndata: Access and store annotated data matrices. Journal of Open Source Software 9(101):4371. DOI 10.21105/joss.04371
- Wolf FA, Angerer P, Theis FJ (2018) SCANPY: large-scale single-cell gene expression data analysis. Genome Biology 19:15. DOI 10.1186/s13059-017-1382-0
- Hao Y, et al. (2024) Dictionary learning for integrative, multimodal and scalable single-cell analysis (Seurat v5). Nature Biotechnology 42(2):293-304. DOI 10.1038/s41587-023-01767-y
- Amezquita RA, Lun ATL, Becht E, et al. (2020) Orchestrating single-cell analysis with Bioconductor. Nature Methods 17(2):137-145. DOI 10.1038/s41592-019-0654-x
- Bredikhin D, Kats I, Stegle O (2022) MUON: multimodal omics analysis framework. Genome Biology 23:42. DOI 10.1186/s13059-021-02577-8
- Lun ATL, Riesenfeld S, Andrews T, et al. (2019) EmptyDrops: distinguishing cells from empty droplets. Genome Biology 20:63. DOI 10.1186/s13059-019-1662-y
