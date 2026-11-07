# Discover or retrieve links

```bash
python scripts/link.py --email researcher@example.org --source gene --target protein --id 672 --linkname gene_protein_refseq --out links.xml
```

The script first prints ESummary labels for source UIDs; confirm they are the intended records. Repeat `--id` for per-source LinkSets. For discovery use `--cmd acheck` and omit target and linkname. Empty linksets may be valid; wrong namespaces can silently resolve to another record. `references/link_catalog.md` lists source-tested linknames; check the live `acheck` output for the chosen source. ELink returns asserted relationships, not sequence similarity.

The source-checking helper is `scripts/check_source_ids.py`; the command above supplies the requester's email to it.

Done when source identities, linkname, direction and returned IDs are checked and the XML saved.
