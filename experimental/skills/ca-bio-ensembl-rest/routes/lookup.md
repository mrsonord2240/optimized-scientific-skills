# Gene, sequence or region lookup

```bash
python scripts/query.py /lookup/symbol/human/BRCA1 gene.json
```

Replace the species and symbol with the actual input; persist the returned stable ID. For expanded gene records, sequence, interval overlap, regulatory features, species and cross-reference endpoints use `references/endpoints.md`. Set `--base` to a validated archive or assembly host when reproducibility requires it. For bulk tables use `bio-biomart-queries`. `scripts/ensembl_client.py` contains reusable retry and domain functions.

Done when the response and its host, release, species and assembly are saved and checked.
