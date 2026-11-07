# Search one database

```bash
python scripts/search.py search --email researcher@example.org --db pubmed --term 'CRISPR[Title]' --retmax 100 --out search.json
```

Choose database, field-qualified term and result limit. The output retains Count, IdList and QueryTranslation; compare Count with returned IDs and do not assume a partial list is complete. Add `--history` to request WebEnv and QueryKey for downstream retrieval; use `bio-batch-downloads` for large downloads. Query design and source-reported pagination limits are in `references/query-design-and-failures.md` and `references/history.md`.

Done when the count, translation, returned IDs or history identifiers and query date are recorded.
