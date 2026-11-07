# Supplementary data or linked runs

```bash
Rscript scripts/supplementary.R GSE123456 geo_output
```

The GEOquery command downloads supplementary files into a fresh directory and writes a manifest; verify completeness. For linked sequencing runs use `python scripts/geo_lookup.py sra --gse GSE123456 --out runs.tsv`; the source pysradb pattern selects the first SRP, so check for additional studies before claiming completeness. Raw run download belongs to `bio-sra-data`. See `references/entrez-search-and-links.md` and `references/geoparse-geoquery.md`.

Done when the requested files or run list are saved, completeness checked and source and processing choice reported.
