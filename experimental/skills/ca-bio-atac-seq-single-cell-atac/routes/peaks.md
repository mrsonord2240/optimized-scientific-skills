# Doublets and cluster peaks

```bash
Rscript scripts/cluster_peaks.R scatac_signac.rds /path/to/macs3 cluster_peaks.rds
```

The input is a clustered Signac object with fragment data; the command calls MACS3 per cluster with the source BED shift/extension settings and saves peaks as RDS. Run doublet detection separately: choose AMULET, ArchR or scDblFinder by depth and input, and check barcode conventions in `references/ecosystem-workflows.md`. Low-depth clusters may need pooling. Use `bio-atac-seq-consensus-peakset` for consensus construction and `bio-atac-seq-differential-accessibility` for replicate-aware differential analysis. Finalize peaks before motif scoring.

Done when doublet handling, per-cluster peak depth and consensus strategy are documented.
