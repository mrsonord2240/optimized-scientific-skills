# Allele-pattern attribution

```bash
python scripts/deconvolute_bystander.py Alleles_frequency_table.zip --target-pos 67 --bystander-pos 69 --out allele_patterns.tsv
```

Positions are 1-indexed within the aligned allele string. Read fractions by pattern alone do not establish which edit caused a phenotype. See `references/screen-analysis.md`.

Done when allele patterns and attribution uncertainty are reported.
