# Subset samples
```bash
bcftools view -S samples.txt input.vcf.gz -Ou | bcftools +fill-tags -Oz -o subset.vcf.gz -- -t AC,AN,AF
bcftools index subset.vcf.gz
```
**Done when:** output samples and recomputed counts match the selected list and output is indexed. One sample name per line; verify plugin and input header. `-s` accepts a comma-separated list, and `-S ^samples.txt` excludes a list. For regions, `bcftools view -r chr1:1000000-2000000 input.vcf.gz -Oz -o region.vcf.gz` requires an index; `-R` uses region files, while `-t/-T` streams targets and uses different overlap semantics. Merge overlapping intervals if version behavior is unknown.
For splitting by sample or chromosome, iterate `bcftools view -s` or `-r` with unique output names and index each output. Rename samples with `bcftools reheader -s rename.txt input.vcf.gz -o renamed.vcf.gz` (old and new names separated by a tab); verify with `bcftools query -l`. Reheader does not rewrite record CHROM: use `bcftools annotate --rename-chrs map.txt input.vcf.gz -Oz -o renamed-contigs.vcf.gz` for a validated old/new contig map. Neither changes reference build. See [details](../references/operations.md).
