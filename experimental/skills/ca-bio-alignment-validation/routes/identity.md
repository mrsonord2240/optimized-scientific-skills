# Contamination and sample identity

```bash
verifybamid2 --SVDPrefix /resources/1000g.phase3.10k.b38.vcf.gz.dat --Reference ref.fa --BamFile sample.bam --Output sample.contam
```

Choose a marker panel compatible with the reference; genome/exome-scale coverage is needed. For relatedness, sex checks, and tumor-normal identity use `references/contamination-swap.md`.

Done when the requested contamination or identity result is interpreted.
