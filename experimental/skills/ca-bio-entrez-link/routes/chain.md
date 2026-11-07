# Chain or score links

```bash
python scripts/link.py --email researcher@example.org --source pubmed --target pubmed --id 35412348 --linkname pubmed_pubmed --cmd neighbor_score --out scored.xml
```

Relevance scores are raw values, not percentages. For gene-to-protein-to-structure or BioProject-to-SRA patterns use `references/code_patterns.md`; resolve accessions to the correct source UID first. Reverse and forward tables can differ. Use `bio-entrez-fetch` to retrieve linked records after confirming their database.

Done when each step's source and target IDs, linkname and direction are recorded with any scores.
