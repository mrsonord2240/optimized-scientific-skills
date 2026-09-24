## Unpaired / Diagonal: GLUE (Python)

**Goal:** Align independent scRNA and scATAC with no shared cells via a prior feature graph.

**Approach:** Configure each dataset with a count-appropriate probabilistic model, build a gene-anchored guidance graph, fit GLUE, then read aligned embeddings.

```python
import scglue

scglue.models.configure_dataset(rna, 'NB', use_highly_variable=True, use_rep='X_pca')     # NB needs RAW counts
scglue.models.configure_dataset(atac, 'ZINB', use_highly_variable=True, use_rep='X_lsi')
graph = scglue.genomics.rna_anchored_guidance_graph(rna, atac)     # peak-near-gene prior; coords must share genome build
# GLUE also trains a VAE, so the seed is pinned explicitly for reproducibility (scglue's
# documented default is already random_seed=0; the totalVI-sized drift was not measured for
# GLUE). Checked against scglue's documented API, not run -- scglue has no Windows build.
glue = scglue.models.fit_SCGLUE({'rna': rna, 'atac': atac}, graph, init_kws={'random_seed': 0})
rna.obsm['X_glue'] = glue.encode_data('rna', rna)
atac.obsm['X_glue'] = glue.encode_data('atac', atac)
```

Verify cell-type structure is preserved (not just modality overlap); adversarial alignment can over-mix distinct populations.

## Unpaired / Diagonal: Seurat v5 Bridge Integration (R)

**Goal:** Map an unpaired scATAC query onto a labelled scRNA reference, using a paired multiome dataset as the bridge.

**Approach:** Preprocess the three datasets in their native pipelines (query ATAC only TF-IDF), build the bridge reference from the scRNA reference plus the multiome bridge, find anchors by projecting the query into the bridge's ATAC LSI space, then transfer labels and project onto the reference UMAP.

```bash
# When Seurat/Signac are in R's normal library paths:
Rscript scripts/seurat_bridge_integration.R rna.rds multi.rds atac.rds bridge_query.rds   # optional: ndims first_lsi_dim SCT|LogNormalize

# Otherwise, expose the private library before starting R (Git Bash/Linux/macOS):
BIO_SKILLS_R_LIB=/path/to/R-lib Rscript scripts/seurat_bridge_integration.R rna.rds multi.rds atac.rds bridge_query.rds
```

In PowerShell, use `$env:BIO_SKILLS_R_LIB = 'C:\path\to\R-lib'` before the `Rscript` command. The script
adds this directory to `.libPaths()` before checking and loading Seurat/Signac, and stops with a direct
setup error if either package is still unavailable.

Inputs, all Seurat objects saved with `saveRDS`: `rna.rds` is the labelled scRNA reference (`meta.data$celltype`; NormalizeData/ScaleData/RunPCA, and `RunUMAP(return.model = TRUE)`); `multi.rds` is the paired multiome bridge with `RNA` (normalized) and `ATAC` (RunTFIDF, RunSVD -> `lsi`) assays; `atac.rds` is the unpaired scATAC query on the SAME peak set as the bridge's ATAC assay, RunTFIDF only. The output query carries `predicted.celltype`, `predicted.celltype.score` and `ref.umap`.

The script starts the ATAC dims at 2 (drops LSI_1): keep that only if DepthCor confirms it tracks depth (see `references/multiome-mofa.md`), otherwise pass `1` as `first_lsi_dim`. It uses `normalization.method = 'LogNormalize'`; pass `SCT` if the reference and bridge RNA were SCTransformed. `MapQuery` must take the object `PrepareBridgeReference` RETURNED as `reference`, not the original scRNA object (the anchorset lives in its `Bridge` assay; passing the original errors "assay ... does not match").

Checked on Seurat 5.5.0 / Signac 1.17.1 with synthetic three-population data (300 reference RNA cells, 300 bridge multiome cells, 300 ATAC-only query cells, planted marker genes and peaks): 81% of query cells recovered their true type with `LogNormalize` (chance is 33%), 75% with `SCT`. Check `predicted.celltype.score` before trusting the transferred labels.
