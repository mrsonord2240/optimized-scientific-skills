# rMATS event plots

```bash
bash scripts/rmats2sashimiplot_events.sh SE.MATS.JC.txt SE sashimi_rmats ctrl1.bam,ctrl2.bam trt1.bam,trt2.bam Control Treatment
```

The wrapper filters FDR < 0.05 and absolute delta PSI > 0.1 by default and checks expected PDFs. Event types and overrides are in `references/rmats2sashimiplot.md`.

Done when all requested event PDFs are present and visually inspected.
