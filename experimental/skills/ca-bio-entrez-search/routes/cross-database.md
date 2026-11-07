# Counts across databases

```bash
python scripts/search.py counts --email researcher@example.org --term 'QUERY' --out counts.json
```

The source default searches ten named databases sequentially, not all Entrez databases. Override using comma-separated `--dbs`; each database requires a request. These are counts, not fetched records. Inspect the per-database query translations; empty results can reflect index lag or differing field semantics.

Done when counts, covered databases, query translations and retrieval date are recorded.
