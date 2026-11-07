# Database fields or spell-check

```bash
python scripts/search.py info --email researcher@example.org --db pubmed --out database_info.json
```

EInfo reports fields, links and update information; omit `--db` to list available databases. For a spelling suggestion use `python scripts/search.py spell --email researcher@example.org --db pubmed --term 'QUERY' --out spelling.json`. Confirm suggestions before changing the user's query. See `references/query-design-and-failures.md` for field and taxonomy-expansion caveats.

Done when requested field metadata or spelling suggestion is saved and interpreted.
