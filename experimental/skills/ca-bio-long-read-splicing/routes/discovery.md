# Isoform discovery and quantification

```bash
isoquant --reference reference.fa --genedb annotation.gtf --fastq sample.fastq.gz --data_type pacbio_ccs --output isoquant_output --threads 16
```

Select data_type and strandedness for the platform; omit genedb for de novo discovery. For ambiguity, novel junction caveats and output names see `references/isoquant.md`. For joint annotation-aware discovery use Bambu and `references/bambu.md`.

Done when requested models and count tables are written and reviewed.
