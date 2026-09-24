## CITE-seq: Denoise ADT, Then Joint Embed (Seurat)

**Goal:** Remove ADT background with DSB before WNN, because WNN does not denoise protein.

**Approach:** Estimate ambient from empty droplets and per-cell technical noise from a mixture plus isotype controls, then feed denoised ADT into the standard PCA -> WNN flow.

Run `examples/cite_seq_analysis.R` for this whole flow (it reads `raw_feature_bc_matrix/` and `filtered_feature_bc_matrix/`, builds `adt_cells` and `adt_empty`, guards the empty-droplet matrix, calls `DSBNormalizeProtein`, then does the PCA -> WNN steps below). Notes on the DSB call:

- DSB gives NO error or warning when a filtered/cell matrix is passed as `empty_drop_matrix` (verified, dsb 2.0.1; see Common Errors). True empty droplets carry mostly ambient signal, so their median total ADT count is markedly lower than in called cells; the example stops if the empty-droplet median is at least 0.5x the cell median.
- `isotype.control.name.vec` must name the ACTUAL isotype rows (often IgG1/IgG2a/Mouse-IgG2b-Ctrl); the example's regex `'[Ii]sotype|IgG'` misses those.
- When isotypes are absent or not matched, set `use.isotype.control = FALSE` (keep `denoise.counts = TRUE`) and pass real names explicitly.

## CITE-seq: WNN Joint Clustering (Seurat)

**Goal:** Build one weighted-NN graph from denoised RNA and ADT and cluster on it.

**Approach:** Reduce each modality independently (PCA on RNA, PCA on the small ADT panel), then learn per-cell modality weights and cluster/embed on the joint graph.

```r
obj[['ADT']] <- CreateAssay5Object(data = adt_dsb)        # DSB output is already normalized data
DefaultAssay(obj) <- 'RNA'
obj <- NormalizeData(obj) |> FindVariableFeatures() |> ScaleData() |> RunPCA(reduction.name = 'pca')

DefaultAssay(obj) <- 'ADT'
VariableFeatures(obj) <- rownames(obj[['ADT']])
obj <- ScaleData(obj) |> RunPCA(reduction.name = 'apca', npcs = min(18, nrow(obj[['ADT']]) - 1))

# dims.list matched to informative dims; small ADT panels saturate by ~1:18
obj <- FindMultiModalNeighbors(obj, reduction.list = list('pca', 'apca'), dims.list = list(1:30, 1:18))
obj <- FindClusters(obj, graph.name = 'wsnn', algorithm = 3)   # algorithm 3 = SLM (the tutorial choice), NOT Leiden
obj <- RunUMAP(obj, nn.name = 'weighted.nn', reduction.name = 'wnn.umap')

# Inspect the per-cell weight distribution; a single dominant modality is a red flag
VlnPlot(obj, features = 'RNA.weight', group.by = 'seurat_clusters')
```
