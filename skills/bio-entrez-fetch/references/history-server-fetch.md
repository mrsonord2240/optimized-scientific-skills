# History-server fetch after ESearch

Read when pulling a large result set (thousands of records) via `webenv`/`query_key`.

Assumes the Required Setup block in `SKILL.md` (`Entrez.email`, imports, `expect_start`).

### History-server fetch (post-ESearch)

**Goal:** Pull a 50,000-record result set without re-sending UIDs.

**Approach:** ESearch with `usehistory='y'`; iterate EFetch with `webenv`/`query_key` and `retstart`. See `batch-downloads` for the production pattern.

```bash
python scripts/history_fetch.py --email you@inst.edu \n  --term 'Homo sapiens[ORGN] AND srcdb_refseq[PROP] AND biomol_mrna[PROP]' --out out.fasta
```

`scripts/history_fetch.py` runs ESearch with `usehistory='y'`, then EFetch in 500-record chunks by `retstart` with `webenv`/`query_key`, sleeping 0.34 s (0.1 s with `NCBI_API_KEY` set). Import `history_fetch(term, out_path, db, chunk)` to call it from Python.
