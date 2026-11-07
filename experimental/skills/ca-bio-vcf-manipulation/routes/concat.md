# Stitch regions and sort
```bash
bcftools concat -f ordered_files.txt -Oz -o genome.vcf.gz
bcftools index genome.vcf.gz
```
**Done when:** output sample order matches all inputs, regions are in reference order, output is sorted and indexed. Use the same samples and order, sorted files and an ordered input list. For overlapping regions, `-a` requires indexed inputs; choose `-d exact` only when intended to remove duplicates. `--naive` requires identical headers and sample order; do not use `--naive-force` without verified identical headers. To sort, run `bcftools sort -T /scratch/tmp -m 4G input.vcf.gz -Oz -o sorted.vcf.gz` before indexing. See [combination details](../references/operations.md).
