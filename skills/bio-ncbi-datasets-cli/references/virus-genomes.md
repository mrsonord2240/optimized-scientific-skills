# Virus genomes

Read when the request is for virus assemblies, virus metadata or virus proteins.

### Virus genomes

**Reference (NCBI Datasets CLI 18.37.0, checked 2026-09-21):**
```bash
# Metadata first (RefSeq only), then download; virus uses `genome taxon`, not `accession`
datasets summary virus genome taxon "Zika virus" --refseq --as-json-lines   | dataformat tsv virus-genome --fields accession,virus-name,length,host-name,release-date

datasets download virus genome taxon "Zika virus" --refseq     --include genome,cds,protein --filename zika.zip --no-progressbar
unzip -q zika.zip -d zika/    # ncbi_dataset/data/{genomic.fna,cds.fna,protein.faa,data_report.jsonl}
```

Default package is `genomic.fna` + `data_report.jsonl`; `--include` adds `cds`, `protein` (and
`annotation`, which yields `annotation_report.jsonl`). Filters: `--refseq`, `--complete-only`,
`--host`, `--geo-location`, `--released-after`, `--lineage` (SARS-CoV-2 only). Field names come from
`dataformat tsv virus-genome --help` (`--fields` accepts quoted `*` wildcards). Live run on
2026-09-21: 2 RefSeq Zika genomes (NC_012532.1, NC_035889.1).
