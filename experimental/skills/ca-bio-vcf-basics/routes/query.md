# Query fields
```bash
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' input.vcf.gz > genotypes.tsv
```
**Done when:** output columns, sample names, allele indices, missing genotypes and phase separators match the input header and records.
Brackets repeat the expression for each sample. Inspect header declarations before selecting INFO or FORMAT fields; multiallelic and per-genotype values are arrays whose declared Number matters. See [field semantics](../references/fields.md).
