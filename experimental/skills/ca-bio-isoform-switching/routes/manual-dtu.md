# Stage-wise differential transcript usage

```bash
Rscript scripts/manual_dtu.R salmon_quant metadata.tsv tx2gene.tsv results.csv
```

Metadata TSV requires unique sample_id and condition; tx2gene TSV requires unique tx and gene columns. Quantifications are salmon_quant/<sample_id>/quant.sf. The source workflow filters with DRIMSeq, tests with DEXSeq, and applies stageR; see `references/manual-dtu.md`. Two per group is exploratory. The alternative satuRn stage-wise implementation is not supplied; see `references/selection.md`.

Done when gene and transcript adjusted results are saved and interpreted.
