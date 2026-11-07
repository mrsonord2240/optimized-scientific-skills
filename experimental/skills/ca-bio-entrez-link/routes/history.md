# Large linked set

```bash
python scripts/history_link.py --email researcher@example.org --source gene --target protein --ids-file gene_uids.txt --linkname gene_protein_refseq --out history.json
```

Use one UID per line. Source code posts chunks of 200, unions the distinct QueryKeys, then requests `neighbor_history`; linking only the final chunk loses earlier IDs. The resulting history set is a union, not a per-input mapping. Validate source namespaces and compare counts before using WebEnv and QueryKey for retrieval. For bounded source-to-record fetch use `bio-batch-downloads`. History lifetime and service limits require current confirmation.

Done when history identifiers and source and linked counts are checked and the union semantics recorded.
