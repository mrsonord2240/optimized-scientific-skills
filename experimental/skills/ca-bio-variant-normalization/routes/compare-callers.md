# Compare callers or cohorts

```bash
bash scripts/compare_callers.sh reference.fa cmp_out gatk.vcf.gz freebayes.vcf.gz
```

Normalizes each input with `scripts/normalize_vcf.sh` (same flags for all), then runs `bcftools isec -p cmp_out/isec`. For two inputs: `0000.vcf` private to the first, `0001.vcf` private to the second, `0002.vcf` and `0003.vcf` shared (records from each file). `ATOMIZE=0` applies to all inputs.

- Every input must be called against the same reference build.
- Normalizing only one side, or using different tools or flag sets per cohort, makes every MNP a false private variant.
- Merge or isec on non-normalized input silently drops matches.

Done when `cmp_out/isec` holds the private and shared sets and each input passed the REF check.
