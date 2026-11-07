# How much needs normalizing

```bash
python scripts/check_normalization.py input.vcf.gz
```

Prints total variants, multiallelic sites and MNPs with cyvcf2. It cannot see indels needing left-alignment (that needs reference context), so the count is a lower bound; normalize with `routes/normalize.md` regardless.

Done when the counts are printed.
