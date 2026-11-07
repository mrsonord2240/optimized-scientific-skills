# Train chromBPNet

```bash
bash scripts/chrombpnet_pipeline.sh atac.dedup.bam peaks.narrowPeak genome.fa chrom.sizes output ATAC 0.5
```

Input peaks must be 10-column narrowPeak. The wrapper makes chromosome splits and non-peaks, trains bias and accessibility models, and runs QC. Choose held-out chromosomes appropriate to the genome using TEST_CHROMS and VALID_CHROMS; read `references/method-reference.md` for inputs, runtime and gate thresholds. Use a fresh output or documented RESUME behavior. Full-scale training was not run in the source evaluation.

Done when the model and held-out QC exist and bias-model motifs have been reviewed.
