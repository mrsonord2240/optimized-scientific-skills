# Search GEO or a publication

```bash
python scripts/geo_lookup.py search --email researcher@example.org --term 'QUERY AND gse[Entry Type]' --out studies.json
```

The query searches Entrez `gds`, which indexes multiple GEO record types. Inspect Count and returned IDs; add organism and technology fields when appropriate. For a publication use `python scripts/geo_lookup.py pubmed --email researcher@example.org --pmid 12345678 --out studies.json`, replacing the example with the requested identifier. Check each returned GSE relation before treating it as a single experiment. Search and link details are in `references/entrez-search-and-links.md`.

Done when matching records, query and limits are saved and each selected Series is inspected.
