# Resolve run accessions

```bash
python scripts/pysradb_resolve.py PROJECT_ACCESSION > runs.txt
```

Replace the placeholder with a project, study, experiment or GSE accession. The source helper reports counts on stderr and run IDs on stdout; verify the full list and metadata before transfer. `bio-entrez-search` and `bio-geo-data` cover related discovery. Source uses pysradb 2.2+ and reported checking 2.5.1; check current version and API.

Done when run IDs are resolved and their access status, expected size and provenance are checked.
