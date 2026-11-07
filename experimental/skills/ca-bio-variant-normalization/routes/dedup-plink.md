# Duplicates and PLINK-ready VCF

```bash
bcftools norm -d exact input.vcf.gz -Oz -o deduped.vcf.gz
```

`-d` modes: `exact` same CHROM, POS, REF, ALT; `snps`; `indels`; `both`; `all` every duplicate at a position; `none` keep (default).

PLINK needs biallelic, SNP-only, deduplicated input:

```bash
bcftools norm -f reference.fa -m- input.vcf.gz | \
    bcftools view -v snps | \
    bcftools norm -d exact -Oz -o gwas_ready.vcf.gz
bcftools index gwas_ready.vcf.gz
```

Done when the output is indexed and no duplicate records remain.
