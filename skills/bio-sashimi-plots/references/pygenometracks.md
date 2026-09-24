## pyGenomeTracks for Multi-Track Figures

**Goal:** Combine splicing with chromatin or coverage tracks for publication figures.

**Approach:** Build coverage bedGraphs and a junction BEDPE from the BAMs, define tracks in an INI file (genes, bedGraph/BigWig, BED, links), then run `pyGenomeTracks --tracks tracks.ini --region ... -o figure.pdf`. pyGenomeTracks 3.9 cannot draw a BAM (`InputError ... can not identify file type`).

```bash
# writes ctrl.bedgraph, trt.bedgraph and junctions.bedpe into the output directory ("." here)
scripts/pgt_tracks.sh . ctrl1.bam,ctrl2.bam,ctrl3.bam trt1.bam,trt2.bam,trt3.bam
```

The script merges one BAM per group and runs `bedtools genomecov -split -bga` on each; `-split` is essential, without it introns are filled with coverage. The `junctions.bedpe` step is described below.

```ini
[gene_models]
file = annotation.gtf
height = 3
title = GENCODE v45
fontsize = 10
file_type = gtf

[ctrl_coverage]
file = ctrl.bedgraph
title = Control
color = #1f77b4
height = 3
min_value = 0
max_value = 200
file_type = bedgraph

[trt_coverage]
file = trt.bedgraph
title = Treatment
color = #ff7f0e
height = 3
min_value = 0
max_value = 200
file_type = bedgraph

[junctions]
file = junctions.bedpe
title = Junctions
height = 5
file_type = links
links_type = arcs
```

Arc height grows with the junction's span, so a short `[junctions]` track crops the widest arc (`height = 2` did on a 3-exon locus whose skipping junction spans 60% of the window; 5 shows it whole): raise `height`, or narrow `--region`, until the widest arc is complete, and look at the figure.

Tracks are scaled independently: set the same `min_value`/`max_value` on both coverage tracks (pick `max_value` from the data) or the two groups are not comparable. A BigWig made from the same `-split` bedGraph works too (`file_type = bigwig`).

The `junctions.bedpe` file must be in **BEDPE format** (6 columns: chr1 start1 end1 chr2 start2 end2 [+ optional score]). Convert from regtools .bed12 junctions (the score is the read count summed over the merged BAM; `-s XS` needs XS-tagged BAMs, otherwise the strand is `?`). `pgt_tracks.sh` does this conversion: it merges every BAM into `all_merged.bam` and indexes it (regtools needs an indexed BAM; on an unindexed BAM it wrote nothing), runs `regtools junctions extract -s XS`, and converts BED12 to BEDPE. In the regtools BED12, column 11 is blockSizes (anchor_left, anchor_right) and column 12 is blockStarts (0, intron_length + anchor_left), so intron start = `$2 + a[1]` and intron end = `$2 + b[2]`; the BEDPE rows are `chr s s+1 chr e-1 e score`.

```bash
pyGenomeTracks --tracks tracks.ini --region chr17:43094000-43125000 -o figure.pdf
```
