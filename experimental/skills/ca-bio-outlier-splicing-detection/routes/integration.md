# Variant and RNA outlier integration

```bash
bcftools query -f '%CHROM\t%POS\t%INFO/SpliceAI\n' spliceai_annotated.vcf > spliceai_raw.tsv
```

Join the same patient's FRASER outliers to predicted variants with `references/variant-integration.md`. Check contig naming; the 1 kb window is a heuristic and delta_max >= 0.2 is the source example. Treat concordance as research prioritization.

Done when candidate matches and filtering criteria are reported.
