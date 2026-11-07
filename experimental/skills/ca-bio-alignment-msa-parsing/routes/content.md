# Alignment content

```bash
python scripts/content.py alignment.fasta results.json
```

The output contains IDs, per-column composition and gap counts, gap counts per row, consensus and conserved positions. Conservation counts gaps in the denominator; consensus defaults to 0.5 and conserved columns to 0.8. Set `--consensus-threshold` or `--conserved-threshold` when requested. See `references/consensus.md` for weighted consensus and ambiguity handling.

For sequence extraction by ID or range, and selecting columns occupied by a reference sequence, read `references/extraction.md`.

Done when the requested content has been extracted from the results file.
