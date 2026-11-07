# Interpret a record
```bash
bcftools view -H input.vcf.gz | head
```
**Done when:** each selected record is interpreted using header declarations, sample ploidy and reference build. Consult [field semantics](../references/fields.md) for QUAL/GQ/PL, AD/DP, cardinality, phase, multiallelic and symbolic alleles, END and gVCF blocks. For cyvcf2 iteration, genotypes and writing, see [Python recipes](../references/python.md). Related shelf skills: `bio-vcf-manipulation`, `bio-variant-normalization` and `bio-vcf-statistics`.
