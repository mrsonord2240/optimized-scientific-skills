# Inspect a GEO Series

```bash
python scripts/geo_series.py relation GSE123456 > relation.txt
```

Replace the accession with the requested GSE. If it is a SuperSeries, inspect its SubSeries and platform composition. Then run `python scripts/geo_series.py matrix GSE123456 --outdir output_dir` with an existing fresh output directory to download the standard matrix and print its processing notes. Missing processing notes mean unknown, not unprocessed. Multi-platform studies may have separate matrices; the script's single standard matrix URL does not enumerate them. See `references/geoparse-geoquery.md` and `references/formats.md`; old ArrayExpress and BioStudies paths are in `references/legacy-and-formats.md`.

Done when Series relationships, platform and processing provenance are reported with the retrieved data.
