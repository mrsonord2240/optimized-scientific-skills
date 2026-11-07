# CITE-seq in Seurat: DSB then WNN

```bash
Rscript scripts/cite_seq_wnn.R raw_feature_bc_matrix/ filtered_feature_bc_matrix/ cite_out isotypes=IgG1,IgG2a
```

Inputs are the 10x raw (unfiltered, with empty droplets) and filtered directories. The script runs RNA QC and PCA, takes empty droplets from the raw matrix, runs `DSBNormalizeProtein`, builds the ADT PCA, runs `FindMultiModalNeighbors`, clusters and runs the UMAP. It writes `cite_out.rds`, `cite_out_weights.csv` (per-cell RNA and ADT weights and cluster) and `cite_out_umap.pdf`. Options: `resolution=0.5`, `adt_dims=18`, `rna_dims=30`.

- DSB gives no error when a filtered matrix is passed as the empty-droplet matrix. The script stops if the empty-droplet median total ADT count is 0.5x the cell median or more.
- `isotypes=` must name the actual isotype rows (for example IgG1, IgG2a, Mouse-IgG2b-Ctrl). The default pattern `isotype|IgG` misses other names; with no match DSB runs without isotype controls.
- WNN does not denoise protein; DSB output is already normalized data, so it is stored as the ADT `data` layer.
- Small ADT panels saturate by about 18 dimensions.
- Read `cite_out_weights.csv`: a single dominant modality is a red flag. Rerun with ADT down-weighted or denoised harder and check the clustering holds.

Done when `cite_out.rds` exists and the answer reports the weight distribution.
