# Convert and index
```bash
bcftools view -Ob -o output.bcf input.vcf.gz
bcftools index output.bcf
```
**Done when:** the new BCF and its CSI index exist and their header and sample list match the source. Use a distinct output path and coordinate-sorted input. For compressed VCF use `bcftools view -Oz -o output.vcf.gz input.bcf` then `bcftools index output.vcf.gz`; CSI is the default, while `bcftools index -t output.vcf.gz` requests TBI. `bgzip -c input.vcf > output.vcf.gz` preserves the plain VCF. Indexing does not establish reference-build compatibility.
