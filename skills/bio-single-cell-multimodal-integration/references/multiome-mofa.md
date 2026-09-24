## Multiome (RNA + ATAC, same cell): Native Pipelines, Then Join

**Goal:** Process each modality in its own statistics before joining, because RNA and ATAC have incompatible distributions.

**Approach:** PCA on RNA, TF-IDF + LSI on ATAC (drop depth-correlated components), then WNN. See scatac-analysis for ATAC QC and the binarization/depth-component caveats.

```r
library(Signac)
DefaultAssay(obj) <- 'RNA'
obj <- NormalizeData(obj) |> FindVariableFeatures() |> ScaleData() |> RunPCA()

DefaultAssay(obj) <- 'ATAC'
obj <- RunTFIDF(obj) |> FindTopFeatures(min.cutoff = 'q0') |> RunSVD()
DepthCor(obj)                                          # diagnose which LSI components track depth

# dims = 2:30 drops LSI_1 ONLY if DepthCor confirms it tracks depth (usually true, not guaranteed)
obj <- FindMultiModalNeighbors(obj, reduction.list = list('pca', 'lsi'), dims.list = list(1:30, 2:30))
obj <- RunUMAP(obj, nn.name = 'weighted.nn', reduction.name = 'wnn.umap')
obj <- FindClusters(obj, graph.name = 'wsnn', algorithm = 3)
```

Merging multiome datasets requires a common peak set: re-quantify all cells against unified peaks, or peak-boundary differences manufacture spurious batch structure. The ATAC gene-activity matrix is an approximation, not measured RNA; do not conflate it with the RNA modality.

## MOFA+ (interpretable shared/specific factors)

**Goal:** Decompose modalities into shared latent factors with per-modality variance explained.

**Approach:** Build a MOFA object from per-modality matrices, set likelihoods to match each data type, run, then interpret factor loadings.

```python
import muon as mu

# likelihoods must match data: gaussian for scaled RNA, bernoulli for binarized ATAC, poisson for counts
mu.tl.mofa(mdata, n_factors=15, outfile='mofa_model.hdf5')   # writes mdata.obsm['X_mofa']
```
