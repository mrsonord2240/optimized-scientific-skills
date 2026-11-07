# Ensembl Compara

```bash
python scripts/query.py compara --symbols BRCA1 TP53 --source human --target mouse --out orthologs.csv
```

Replace symbols, species and output path. This uses the source batch client, which retains per-symbol found, no-call and request-failed statuses and sleeps 0.07 seconds between calls. Inspect failed and empty rows; do not treat a partial batch as complete. Filter `type` for `ortholog_one2one` if needed; source records report type, taxonomy level and source/target identity, while confidence may be absent. `examples/compara_orthologs.py` also shows symbol resolution, rename handling and batch use. For larger exports consider `bio-biomart-queries`; general Ensembl work uses `bio-ensembl-rest`.

Done when each input has a reviewed status and requested calls have species, type and identifiers recorded.
