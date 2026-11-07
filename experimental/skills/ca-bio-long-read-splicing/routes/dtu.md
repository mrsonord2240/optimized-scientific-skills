# Differential transcript usage

```bash
Rscript scripts/dtu.R reads_manifest.tsv flair_quantified.counts.tsv results.csv gene_ids_without_underscores
```

The no-header manifest has sample_id, condition, batch and reads columns. The required final argument confirms the source workflow's FLAIR ID assumption: gene IDs have no underscores. Validate the mapping before running. Ambiguous and truncated reads leave assignment uncertainty; see `references/dtu.md`.

Done when stage-wise gene and transcript results are saved and interpreted.
