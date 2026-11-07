# Reading the numbers

No script. Expected values are conventions that shift with capture kit, ancestry, reference build and caller; start from them and confirm against a matched cohort.

| Metric | WGS | WES | A deviation means | Action |
|--------|-----|-----|-------------------|--------|
| Ti/Tv, overall | ~2.0-2.1 | ~3.0-3.3 | Low: false-positive transversions (random errors give ~0.5). High: over-filtering removed transversions | Tighten site filters if low; check the filter for transversion bias if high |
| Ti/Tv, novel only | near overall | near overall | Well below the known-site value: false positives concentrate in the novel fraction | Raise stringency |
| Het/hom ratio | ~1.5-1.6 EUR, ~2.0+ AFR | same | High versus same-ancestry peers: contamination or reference bias. Low: inbreeding, ROH or chromosome loss | Infer ancestry first (peddy or somalier), then flag outliers within each group; never one global cutoff |
| Novel fraction (vs dbSNP) | low for common, high for rare | same | Novel common variants: almost always artifacts | Stratify novel% by allele frequency |
| Call rate per sample | >95-98% | >95-98% | Low coverage or poor sample | Drop worst samples, recompute per-variant missingness, iterate |
| Call rate per variant | >95-99% | >95-99% | Hard-to-genotype site | Drop worst variants after sample QC; rare variants tolerate more missing |
| Het allele balance | centered 0.5 | 0.5 | Shifted off 0.5 with elevated het count: contamination signature | `routes/contamination.md` |
| Excess het / HWE | in equilibrium within ancestry | same | Excess het: collapsed paralog or CNV mapping artifact. Het deficit: often real (substructure, inbreeding) | `routes/missingness-hwe.md` |
| Kinship | matches manifest | matches manifest | Unexpected high kinship: swap or duplicate; a "replicate" that is not related: mislabel | `routes/identity.md` |

A WES callset with Ti/Tv ~2.1 means the filtering is too loose.

## Novel versus known

```bash
bcftools annotate -a dbsnp.vcf.gz -c ID input.vcf.gz -Oz -o annotated.vcf.gz
# novel fraction = records with ID "." over total; stratify by INFO/AF
bcftools view -H annotated.vcf.gz | awk '{n++; if($3==".") novel++} END{print "novel:", novel/n}'
```

A single novel% without frequency stratification is uninformative.

Done when each flagged metric has a stated cause and action, compared against a matched cohort.
