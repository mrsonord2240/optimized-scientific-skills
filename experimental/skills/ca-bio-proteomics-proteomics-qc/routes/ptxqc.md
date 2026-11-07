# PTXQC report from a MaxQuant txt/ folder

```bash
Rscript scripts/ptxqc_report.R <maxquant_txt_dir> files_written.txt
```

`createReport(txt_folder = ...)` takes the `txt/` directory (with `summary.txt`), not a data frame. PTXQC is on CRAN, not Bioconductor. It writes per-metric scores in [0,1], a heatmap PDF, mzQC and YAML into that folder, and an HTML report only when Pandoc is reachable. The Pandoc message is not fatal; for HTML, install Pandoc or set `RSTUDIO_PANDOC` and confirm `rmarkdown::pandoc_available()`.

- Scores are against a representative file and generic marks (% MS2 identified under 20% bad, 20-35% ok, 35% and above great; at least 75-85% of peptides at 0 missed cleavages). They are not biological ceilings.
- The general contaminant score has no 1% cutoff; that belongs to PTXQC's user-defined special-contaminant plot. Judge contaminants against the lab baseline.
- PTXQC does not replace `routes/raw-signal.md`; it does not apply the inspect-before-normalize rule.

Done when the report files exist and the low-scoring metrics are named.
