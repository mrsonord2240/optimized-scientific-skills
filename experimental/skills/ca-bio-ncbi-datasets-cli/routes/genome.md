# Genome data

```bash
datasets download genome accession GCF_000001405.40 --include genome,gff3,gtf,protein,cds,seq-report --filename human_grch38.zip
```

Replace the accession and output name as needed. Unzip the archive and check the expected accession and files. For metadata, run `datasets summary genome taxon "Salmonella enterica" --assembly-level chromosome,complete --released-after 2024-01-01 --as-json-lines | dataformat tsv genome --fields accession,organism-name,assminfo-level,assmstats-scaffold-n50,assminfo-release-date > sal_2024.tsv`. Check `dataformat tsv genome --help` for installed field names. `--reference` restricts to reference assemblies; omit it for the full set. The source Python wrapper is available at `scripts/datasets_wrapper.py`; run `python scripts/datasets_wrapper.py --taxon "Escherichia coli" --accession GCF_000005845.2 --out ecoli.zip --include genome,gff3,protein` or import its functions. See [background](../references/background.md).

Done when the archive or summary contains the requested records and expected file types; zero records can still return exit status zero.
