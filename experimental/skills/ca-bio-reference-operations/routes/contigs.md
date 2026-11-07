# Resolve contig names

```bash
bash scripts/rename_contigs.sh sample.bam map.tsv renamed.bam reference.fa
```

The map has old and new names in tab-separated columns; the optional fourth argument checks names and lengths against the reference FAI. Only rename when sequences match: compare lengths and M5 where available. The script reheaders and indexes BAM; SA, XA and OA tags can retain old contig strings. Identify reference flavour, derive a map and handle these tags using `references/contig-naming.md`.

Done when output names and lengths match the intended reference and sequence identity has been checked.
