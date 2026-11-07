# Bulk chromVAR

```bash
Rscript scripts/chromvar_bulk.R peaks.bed counts.tsv depth.tsv out samples=samples.tsv ref=control
```

- `peaks.bed`: chr, start0, end, same order as the rows of `counts.tsv` (peaks x samples).
- `depth.tsv`: header line, then `sample<TAB>total reads or fragments in the library`, same unit as the counts.
- `samples.tsv` (header, then `sample<TAB>condition`) turns on the limma contrast. `ref=` is required with it. Leave both off for variability only (the right choice for 3-5 samples).
- Options: `genome=BSgenome.Hsapiens.UCSC.hg38` (set the BSgenome matching the peaks), `pfm_rds=<PFMatrixList .rds>` for CIS-BP or plant motifs instead of JASPAR2024 vertebrates CORE, `seed=2024`, `niterations=50`, `alpha=0.05`, `min_logfc=0.5`.

The script adds GC bias, drops samples under 1500 reads or 0.15 in-peaks fraction, filters peaks (non-overlapping, 10 fragments), matches motifs at `p.cutoff=5e-5`, samples background peaks, computes deviations, variability and the limma contrast, and writes `out/`: `chromvar_deviations.csv` (z-scores, motifs labelled `ID (name)`), `chromvar_variability.csv`, `chromvar_differential.csv`, two PDFs and `run.log` (seed, genome, motif source).

- `logFC` is the z-score difference between groups; `|logFC| >= 0.5` is a demonstration filter, apply the request's own cutoff to the table.
- With JASPAR2024 the script reads the SQLite handle (TFBSTools does not dispatch on the JASPAR2024 object directly).
- Use `adj.P.Val < 0.05` for differential motifs. Raw t-tests on z-scores are mis-calibrated, and chromVAR is not a differential peak-count method (use differential accessibility).
- Variability needs no condition labels and suits unsupervised TF discovery.

Done when `out/chromvar_variability.csv` exists (and `chromvar_differential.csv` for a contrast) and the answer reports seed, motif database, and peak and sample counts after filtering.
