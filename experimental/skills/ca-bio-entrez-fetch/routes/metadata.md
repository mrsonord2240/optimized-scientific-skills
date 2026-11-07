# Metadata and literature

```bash
python scripts/fetch.py --email researcher@example.org --db pubmed --id 35412348 --rettype medline --retmode text --out article.medline
```

Use `--summary` for ESummary metadata instead of a full EFetch record. For database-specific formats and parsing read `references/pubmed.md`, `references/gene-taxonomy-gds.md`, or `references/sra.md`; SRA runinfo is metadata, not raw reads. Inspect expected fields and schema. Organism inferred from a title is a heuristic; taxonomy IDs provide a more reliable lookup.

Done when requested records are saved and identifiers and fields are verified.
