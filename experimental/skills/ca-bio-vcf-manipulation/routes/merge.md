# Merge samples
```bash
bcftools merge -l files.txt -Oz -o cohort.vcf.gz
bcftools index cohort.vcf.gz
```
**Done when:** indexed output has expected samples, contigs, representations and genotype missingness. List one sorted, indexed input per line; use distinct output paths. Normalize inputs consistently against the same reference first. Merge unions samples, not regions, and is not joint genotyping: absence of a site is not evidence for hom-reference. Default missing genotypes remain missing; `-0` forces hom-reference and is justified only with confirmed coverage. Ensure sample names and record contigs match; `--force-samples` renames collisions rather than resolving their identity. See [combination details](../references/operations.md).
