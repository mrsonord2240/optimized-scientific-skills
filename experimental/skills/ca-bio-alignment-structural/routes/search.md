# Structural search

```bash
foldseek easy-search query.pdb database result.m8 tmp/
```

Use the local 3Di+AA default for broad search. For `--alignment-type 1`, filter on alignment TM-score rather than E-value; do not treat a row count equal to `--max-seqs` as the total number of hits. See `references/search.md` for custom fields, clustering and low-pLDDT masking. For complex searches use `references/foldseek-multimer.md`; for weak fold matches consider DALI in `references/dalilite.md`.

Done when the hit file is produced and reported with method, filters and coverage.
