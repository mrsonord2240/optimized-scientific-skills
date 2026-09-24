## rmats2sashimiplot

**Goal:** Plot directly from rMATS event coordinates without manual region calculation.

**Approach:** Filter the rMATS event file to the events to plot (it draws every row), pass BAM lists + a group file + event type, then check the output: rmats2sashimiplot **exits 0 when it fails** and leaves `Sashimi_plot/` empty. The contig is matched to the BAM header automatically (`chrX` in the event file, `X` in the BAM works).

```bash
scripts/rmats2sashimiplot_events.sh rmats_output/SE.MATS.JC.txt SE sashimi_rmats \
    ctrl1.bam,ctrl2.bam,ctrl3.bam trt1.bam,trt2.bam,trt3.bam Control Treatment
```

The script keeps only the significant events (`FDR` < 0.05 and |`IncLevelDifference`| > 0.1, columns found by header name; env `FDR`, `DPSI`), writes the group file (`label: first-last`, 1-based over the `--b1` replicates then the `--b2` replicates), runs rmats2sashimiplot with `--exon_s 1 --intron_s 5 --group-info --color '#1f77b4,#ff7f0e'` (env `COLORS`) and exits 1 unless `sashimi_rmats/Sashimi_plot` holds one non-empty PDF per event.

`--event-type` (4.0.0; the old `-t SE` is rejected, rc 2) takes SE, A5SS, A3SS, MXE or RI. `--exon_s 1 --intron_s 5` draws introns at 1/5 of their real length. `--group-info` gives one plot per group (arc labels = group mean, plus the group's mean IncLevel); without it there is one plot per replicate and `--color` needs one colour per replicate, otherwise it prints `Error: Must provide sample label and color for each entry in bam_files!` and still exits 0.
