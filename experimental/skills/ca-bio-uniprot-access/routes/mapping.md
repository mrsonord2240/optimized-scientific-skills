# Identifier mapping

```bash
python scripts/query.py map --ids ID1 ID2 --from-db Ensembl --to-db UniProtKB --out mapping.json
```

Replace IDs and namespaces with the desired input and target. The source client submits an asynchronous job, polls with a 600-second default timeout and three-second interval, reads the full result stream, and reports `failedIds`. Check all failures, multiplicity and reviewed status; mapping can be ambiguous. For obsolete accessions set `--from-db UniProtKB_AC-ID`. Confirm currently supported namespace names in the service catalog. `bio-biomart-queries` is an alternative for large Ensembl-based batches.

Done when all input identifiers have mapped results or explicit failures and the resource release is recorded.
