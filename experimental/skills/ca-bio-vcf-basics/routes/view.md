# Inspect VCF
```bash
bcftools view -h input.vcf.gz
```
**Done when:** the header, declared fields, contigs and samples have been reviewed against the intended reference build.
Use `bcftools view -H input.vcf.gz` for records, or `bcftools query -l input.vcf.gz` for sample names. `bcftools view -r chr1:1000-2000 input.vcf.gz` requires coordinate-sorted bgzip-compressed indexed input and contig names matching its header.
For a cyvcf2 summary run `python scripts/view_vcf.py input.vcf.gz 10`. The optional positional limit defaults to 10. Its FILTER summary labels missing and PASS together; consult the raw record to distinguish them. See [field semantics](../references/fields.md) and [Python recipes](../references/python.md).
