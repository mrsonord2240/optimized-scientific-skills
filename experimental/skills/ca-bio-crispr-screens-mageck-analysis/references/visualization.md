## MAGeCKFlute Integration (R)

**Goal:** Use the canonical pathway-analysis dashboard for downstream visualization.

**Approach:** Load MAGeCK output in R, run FluteRRA or FluteMLE, which produces volcano, rank, square-plot, and KEGG/Reactome enrichment.

```r
library(MAGeCKFlute)
# FluteRRA/FluteMLE do NOT create their own output directory -- outdir must
# already exist, and the real subfolder they write into is outdir/MAGeCKFlute_<proj>/,
# not outdir/ itself. proj defaults to NA, which still needs the matching
# MAGeCKFlute_NA/ directory pre-created -- always pass proj explicitly instead.
dir.create("flute_output/MAGeCKFlute_test1", recursive = TRUE, showWarnings = FALSE)

# After mageck test
FluteRRA(gene_summary = "drug_vs_veh.gene_summary.txt",
         sgrna_summary = "drug_vs_veh.sgrna_summary.txt",
         proj = "test1",
         organism = "hsa",
         outdir = "flute_output/")

# After mageck mle
dir.create("flute_mle_output/MAGeCKFlute_test1", recursive = TRUE, showWarnings = FALSE)
FluteMLE(gene_summary = "timecourse_mle.gene_summary.txt",
         treatname = "day21", ctrlname = "day7",
         proj = "test1",
         organism = "hsa",
         outdir = "flute_mle_output/")
```

For `FluteMLE`, both `treatname` and `ctrlname` must be condition names emitted in the MLE
`gene_summary` (for example `day7` and `day21`). The all-ones `baseline` design-matrix column is the
model intercept; it is not an emitted condition and cannot be used as `ctrlname`.

## MAGeCK-VISPR Interactive Dashboard

**Platform note:** `mageck-vispr` ships only via bioconda for Linux/Mac (Snakemake + Flask); there is
no PyPI or Windows build. On Windows, use bio-crispr-screens-screen-qc's pandas/matplotlib QC instead -- it computes
the same per-sample Gini/mapping-rate/count-summary numbers this dashboard visualizes, without the
interactive server.

```bash
mageck-vispr init my_screen
# Edit config.yaml to point to counts + library
snakemake --cores 8                 # mageck-vispr generates a Snakemake workflow; run it with snakemake
vispr server results/*.vispr.yaml   # serve the interactive dashboard
```
