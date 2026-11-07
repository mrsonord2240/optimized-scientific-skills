# Resumable query retrieval

```bash
python scripts/history.py nucleotide 'QUERY' output.fasta checkpoint.json --email researcher@example.org
```

Use a fresh output path or a checkpoint from the same query and output; the source worker truncates to the recorded byte boundary when resuming. Save the query, database and snapshot context: a changing search result order can make an offset-based resume inconsistent. Set `NCBI_API_KEY` in the environment if available. Source default is FASTA, 500 records per batch and five retries. Verify count and integrity after completion. For rate, retry and history caveats read `references/history.md`. For literal entire databases or over one million sequences, use NCBI bulk mirrors rather than an EFetch loop.

The wrapper calls `scripts/robust_download.py`; do not run that worker directly with its hardcoded demonstration query.

Done when the requested records are downloaded, count and integrity are checked and query and retrieval date recorded.
