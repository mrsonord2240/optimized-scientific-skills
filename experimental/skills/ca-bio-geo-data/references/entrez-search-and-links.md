# GEO search and linked runs

Use `scripts/geo_lookup.py search` with a field-qualified GDS query or `pubmed` to retrieve linked GDS summaries. Inspect result Count and IDs, then check the selected GSE's SOFT relation with `scripts/geo_series.py`.

The `sra` operation uses the source pysradb GSE-to-SRP then SRP-to-SRR pattern. It selects the first returned SRP and prints the complete SRP list; inspect the other studies when present. The source also documented a GDS-to-BioProject-to-SRA ELink fallback, but it is not implemented by this entry point. Use `bio-entrez-link` for deliberate cross-database navigation and `bio-sra-data` for run downloads.
