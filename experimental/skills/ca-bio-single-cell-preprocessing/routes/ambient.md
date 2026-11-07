# Ambient RNA removal

Pick ONE tool; stacking tools compounds over-removal. Ambient removal runs on the raw matrix, before QC.

| Tool | Needs | Command |
|------|-------|---------|
| SoupX (R) | Cell Ranger `outs/` with BOTH raw and filtered matrices | `Rscript scripts/soupx.R cellranger_outs/ adjusted_counts.rds` |
| DecontX (R, `celda`) | Filtered cells only (no raw needed); relies on cluster purity | `celda::decontX(counts)` on a SingleCellExperiment |
| CellBender (Python CLI, normally GPU) | RAW `.h5`; also calls cells; strong on nuclei | `cellbender remove-background --input raw_feature_bc_matrix.h5 --output out.h5 --fpr 0.01` |

- SoupX needs clusters. If Cell Ranger's clustering output is absent the script derives them with a quick Seurat pass and stops if it cannot; do not accept an unclustered estimate. `autoEstCont` fails on homogeneous data.
- The adjusted counts are rounded to integers (`roundToInt=TRUE`) for negative-binomial models.
- After any tool, check that a known-specific marker survives and that off-target markers (Hb, Ig) fell.
- SoupX subtracts a per-cell scalar of one global soup profile; CellBender learns a per-droplet background. They disagree and neither is proven better; CellBender over-removes real low-abundance genes at high `--fpr`. Matters most for solid tumors, snRNA-seq and blood.
- The CellBender and DecontX commands above are not wrapped in a script.

Done when the contamination fraction (SoupX rho) or the tool's equivalent is reported and the adjusted matrix exists.
