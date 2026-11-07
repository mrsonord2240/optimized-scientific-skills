# Doublets and cluster peaks

```bash
Rscript scripts/signac_workflow.R filtered_peak_bc_matrix.h5 fragments.tsv.gz singlecell.csv
```

Run doublet detection separately: choose AMULET, ArchR or scDblFinder by depth and input, and check barcode conventions in `references/ecosystem-workflows.md`. For cluster pseudobulk MACS3 and ArchR peak calls see the same reference; low-depth clusters may need pooling. Use `bio-atac-seq-consensus-peakset` for consensus construction and `bio-atac-seq-differential-accessibility` for replicate-aware differential analysis. Finalize peaks before motif scoring.

Done when doublet handling, per-cluster peak depth and consensus strategy are documented.
