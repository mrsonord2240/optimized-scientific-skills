# Coordinates and release pinning

```bash
python scripts/query.py --dataset hsapiens_gene_ensembl --attributes ensembl_gene_id,external_gene_name,chromosome_name,start_position,end_position,strand,gene_biotype --filter chromosome_name --value 17 --output coordinates.tsv
```

Replace chromosome and dataset with the requested values; `biotype` is a filter and `gene_biotype` an attribute in the source schema. This Python route targets the current Ensembl host; for release pinning use the biomaRt example in `references/query-patterns.md`. Check coordinate convention and genome assembly. Chunk large results by chromosome and inspect for truncation.

Done when the requested table, build and release are verified and saved.
