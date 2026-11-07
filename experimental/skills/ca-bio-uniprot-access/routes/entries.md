# Entry and cross-references

```bash
python scripts/query.py entry --accession ACCESSION --out entry.json
```

Replace the placeholder with an accession. The source client extracts key fields from nested entry JSON and raises for inactive entries. Use `python scripts/query.py xrefs --accession ACCESSION --out xrefs.json` to group references, or `python scripts/query.py isoforms --accession ACCESSION --out isoforms.json` to list alternative products. The source client also exposes `fetch_isoform_fasta` for a suffixed accession; a bare accession refers to canonical sequence. Schema fields can be absent and redirects can resolve secondary IDs.

Done when the record or explicit no-result/error is checked and the release and identifier are recorded.
