# Compare normalized callsets
```bash
bcftools isec -p comparison -Oz caller1.norm.vcf.gz caller2.norm.vcf.gz
```
**Done when:** the four output partitions and their source annotations are reviewed against the indexed, consistently normalized input callsets. Same reference build and compatible allele representation are prerequisites; normalization is covered by `bio-variant-normalization`. Default matching considers alleles; position-only matching can conflate alleles. The shared partitions retain their respective source annotations.
For counts, `bash scripts/compare_vcfs.sh caller1.norm.vcf.gz caller2.norm.vcf.gz comparison` accepts two inputs and an optional directory (default `vcf_comparison`); it creates missing indexes beside inputs and requires bcftools, wc, tr and bc. Use a fresh output directory. For selection flags and position-only Python example see [comparison details](../references/operations.md) and [position example](../references/position.md). Structural variant breakpoint matching and gVCF joint genotyping are upstream-only topics here; do not treat ordinary small-variant isec as either.
