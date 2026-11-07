# Variant records

```bash
python scripts/variant_records.py --email researcher@example.org --db clinvar --uid 4887763
```

Use `--db snp` and the numeric SNP UID for dbSNP. The source parser reports a subset of fields; read `references/variants-clinvar-snp.md` for formats and XML limitations. Report the record's own classification with its accession and context; it is not a diagnosis.

Done when the requested record or stated classification is retrieved and its source and limitations recorded.
