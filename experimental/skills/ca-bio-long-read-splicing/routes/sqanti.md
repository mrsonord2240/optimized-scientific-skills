# Isoform classification and artifact filtering

```bash
sqanti3_qc.py --isoforms isoforms.gtf --refGTF annotation.gtf --refFasta reference.fa -o sqanti3 -d sqanti3_qc --aligner_choice minimap2 --cpus 8
```

Provide transcript/exon-only GTF and inspect intra-priming and RT-switching flags. Filter and optional TSS/polyA inputs are described in `references/sqanti.md`.

Done when classification and artifact decisions are recorded.
