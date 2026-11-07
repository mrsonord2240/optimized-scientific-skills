# Virus genomes

```bash
datasets download virus genome taxon "Zika virus" --refseq --include genome,cds,protein --filename zika.zip --no-progressbar
```

Replace the taxon and output. For metadata use `datasets summary virus genome taxon "Zika virus" --refseq --as-json-lines | dataformat tsv virus-genome --fields accession,virus-name,length,host-name,release-date`. Source filters include `--complete-only`, `--host`, `--geo-location`, `--released-after`, and SARS-CoV-2 `--lineage`; check installed help. See [virus details](../references/virus-genomes.md).

Done when the downloaded package or metadata contains the intended records and files.
