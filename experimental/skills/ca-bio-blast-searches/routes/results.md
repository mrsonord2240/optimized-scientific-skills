# Inspect saved BLAST XML

```bash
python scripts/parse_results.py hits.xml hits.tsv
```

The table reports the first HSP per alignment, score, E-value, identity and the source's alignment-length divided by query-length coverage measure. Gaps and translated searches complicate that measure; inspect query coordinates and units before applying coverage thresholds. E-values depend on database size; record search space and use bit score carefully across scoring settings. See `references/statistics.md` and `references/results-and-rid.md`.

Done when results are filtered using stated identity, coverage and significance criteria and their limitations reported.
