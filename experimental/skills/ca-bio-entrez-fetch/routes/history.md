# Fetch a query result set

```bash
python scripts/history_fetch.py --email researcher@example.org --term 'QUERY' --out output.fasta
```

The source script writes a fresh FASTA in chunks of 500 and does not checkpoint or validate completeness; choose a new output path. For resumable large downloads use `bio-batch-downloads`. See `references/history-server-fetch.md` for the underlying history workflow and rate behavior.

Done when the expected count and output integrity are checked and query and retrieval date recorded.
