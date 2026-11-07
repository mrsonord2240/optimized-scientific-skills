# Multi-track figure

```bash
bash scripts/pgt_tracks.sh . ctrl1.bam,ctrl2.bam trt1.bam,trt2.bam
```

The wrapper generates group bedGraphs and BEDPE junctions; configure tracks and run pyGenomeTracks as in `references/pygenometracks.md`. Use matching coverage scales and enough arc height. BAM is not a direct coverage track in this workflow.

Done when the figure shows complete arcs and comparable coverage scales.
