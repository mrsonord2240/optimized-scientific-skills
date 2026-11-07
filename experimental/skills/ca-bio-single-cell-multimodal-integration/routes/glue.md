# GLUE (unpaired scRNA and scATAC)

```bash
python scripts/glue_unpaired.py rna.h5ad atac.h5ad glue_out
```

`rna.h5ad`: raw counts in `X` (NB likelihood), highly variable genes flagged, `X_pca` in `obsm`, gene coordinates in `var`. `atac.h5ad`: counts (ZINB), `X_lsi`, peak names `chr:start-end`. RNA and ATAC coordinates must share one genome build. Output: `glue_out_rna.h5ad` and `glue_out_atac.h5ad` with `obsm['X_glue']`; the script prints the guidance-graph edge count.

- A near-zero edge count means a genome-build or coordinate mismatch; the alignment is then meaningless. Align coordinates to one build first.
- Check that cell-type structure survives, not only modality overlap: adversarial alignment can over-mix distinct populations.
- `scglue` does not install on Windows (pybedtools needs pysam, no wheel): use WSL, Linux or macOS. The script is checked against the documented API only, not run.
- The seed is pinned (`random_seed: 0`).

Done when both output files exist and the answer reports the edge count and a cell-type preservation check.
