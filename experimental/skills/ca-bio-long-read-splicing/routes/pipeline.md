# Annotated long-read pipeline

```bash
PLATFORM=ont REFERENCE=reference.fa GTF=annotation.gtf FASTQ=reads.fastq.gz SAMPLE=s1 bash scripts/longread_splicing_pipeline.sh
```

PLATFORM is hifi, ont (unstranded cDNA), or drna (direct RNA). Set SR_JUNCTIONS to STAR junctions when available; inspect outputs and SQANTI filters. This workflow runs both IsoQuant and FLAIR and uses an annotation. For tool disagreements and troubleshooting see `references/interpretation.md`.

Done when alignment, isoform models and filtered classification are reviewed.
