# Compara homology

```bash
python scripts/query.py /homology/symbol/human/BRCA1 orthologs.json --param type=orthologues --param target_species=mouse
```

Replace symbol and species. Preserve orthology type and both identities; missing confidence means unknown, not low confidence. For paralogs and gene trees use `references/endpoints.md`. Record the release and target species before comparing outputs.

Done when requested homologies and their type and available confidence are reported.
