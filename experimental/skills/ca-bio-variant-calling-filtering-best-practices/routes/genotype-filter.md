# Genotype-level filters

```bash
bash scripts/genotype_filter.sh site_filtered.vcf.gz gt_filtered.vcf.gz gq=20 dp=8
```

Sets genotypes to `./.` and keeps the site. Defaults `gq=20` (below 99% confidence) and `dp=8` (use 8-10 for WGS); `gq=0` or `dp=0` skips that filter. Add `ab=1` to also null hets whose alt fraction from AD is 0.2 or less, or 0.8 or more (a true het is near 0.5; far from it suggests a mapping artifact, CNV or contamination). GATK does not emit allele balance, so it is derived from AD.

- Run on a site-filtered VCF, and before cohort missingness, HWE or allele frequency. HWE on a matrix full of low-GQ genotypes manufactures deviation.
- Allele balance must be a genotype-level `-S .` filter. A site-level `-i 'GT="het" & ...'` deletes every site with no het genotype, including all hom-alt-only sites.
- The script prints how many genotypes became no-call; zero means the settings changed nothing.
- `./.` is not `0/0`: only no-calls count as missing.

Done when the printed "newly set" count is reported and the record count equals the input.
