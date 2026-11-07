# Search or stream metadata

```bash
python scripts/query.py search --query 'reviewed:true' --fields accession,id --out results.tsv
```

Replace the query and fields with the intended scope. Search returns at most one page, with a warning if total exceeds the page; use `python scripts/query.py stream --query 'reviewed:true' --fields accession,id --out all-results.tsv` only after estimating the size. The client holds the streamed response in memory. UniProt query syntax and field names differ from Entrez; inspect the service field catalog. `reviewed:true` selects reviewed entries; omitting it also includes unreviewed records.

Done when the returned count and requested scope are checked, including any truncation warning.
