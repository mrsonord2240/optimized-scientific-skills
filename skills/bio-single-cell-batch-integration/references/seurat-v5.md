# Seurat v5 integration details

Run [`../scripts/seurat_v5_integration.R`](../scripts/seurat_v5_integration.R) with an `.rds` Seurat object, output path, batch metadata key, and method. The script uses layer splitting, RPCA by default, and sets `future.globals.maxSize` to 4 GiB before anchor integration because Seurat's 500 MiB default fails on routine multi-thousand-cell objects.

`CCAIntegration`, `RPCAIntegration`, and `HarmonyIntegration` are supplied by Seurat. `FastMNNIntegration` and `scVIIntegration` require SeuratWrappers; `scVIIntegration` also requires configured reticulate/scvi-tools. The latter methods should be selected only after those dependencies resolve, not passed as bare symbols to Seurat.
