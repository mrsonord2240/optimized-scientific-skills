# Convert between Python and R

First record what the source holds, so you can diff after:

```bash
python scripts/h5ad_inventory.py data.h5ad
```

Run it again on the converted-and-reloaded file and compare. Cross-ecosystem hops drop `layers`, `obsp`/`varp`, nested `uns`, and turn categoricals into character or NA.

| Direction | Tool | Notes |
|-----------|------|-------|
| h5ad to Seurat or SCE, R only | anndataR | Needs R 4.5 or later; not installable on R 4.4 |
| h5ad to SCE | zellkonverter `readH5AD` | Needs a Python anndata; its default environment can be too old to read files written by anndata 0.13 |
| h5ad to Seurat or SCE, read only | schard `h5ad2sce` | Pure R; reads `X` only, no `counts` layer |
| SCE to h5ad | zellkonverter `writeH5AD` | |
| Live Python plus R session | anndata2ri (rpy2) | |
| Anything to AnnData, one call | sceasy | Set `drop_single_values=FALSE`; `TRUE` deletes every obs/var column with one unique value |
| SeuratDisk | any | Avoid: unmaintained, broken on Seurat v5 |

```r
# h5ad to SCE, X as stored
library(schard)
sce_hvg <- schard::h5ad2sce('data.h5ad', use.raw = FALSE)

# Recover .raw (full gene set). zellkonverter's raw=TRUE can leave altExp empty; schard is the reliable route
sce_raw <- schard::h5ad2sce('data.h5ad', use.raw = TRUE)
length(SingleCellExperiment::altExpNames(sce_hvg))   # 0 means zellkonverter's raw=TRUE silently failed

# Seurat to h5ad (no tool converts this directly)
library(Seurat); library(zellkonverter)
sce <- as.SingleCellExperiment(seurat_obj)   # counts -> counts, data -> logcounts, reductions -> reducedDims
writeH5AD(sce, 'out.h5ad')
```

- `as.SingleCellExperiment()` silently drops `scale.data`. Save it first if needed (it is HVG genes x cells, dense): `saveRDS(as.matrix(LayerData(seurat_obj, layer = 'scale.data')), 'scale_data.rds')`. Reductions such as PCA survive into `obsm`.
- zellkonverter maps `obsm` to `reducedDims`, `varm` to a `rowData` matrix column (not reducedDims), `obsp`/`varp` to `colPairs`/`rowPairs`, and `uns` to `metadata()` (lossy).
- After conversion check the transpose: cells must be columns in R, rows in Python, with `obs` and column metadata swapped accordingly.

Done when the inventory before and after is compared and every slot that disappeared is named.
