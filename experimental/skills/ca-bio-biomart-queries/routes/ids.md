# Bulk ID mapping

```bash
python scripts/query.py --dataset hsapiens_gene_ensembl --attributes ensembl_gene_id,external_gene_name,hgnc_id,refseq_mrna,uniprotswissprot --filter ensembl_gene_id --values-file ids.txt --output mapping.tsv
```

Supply one ID per line. The source helper bypasses pybiomart 0.2.0 ID-list filter discovery, batches from 500 values and halves a batch on HTTP 414. It accepts one large list filter; validate missing IDs and many-to-many rows. Source-tested Ensembl limits certain external-reference attributes per query; split and join if rejected. For GO use `external_gene_name` and the attributes in `references/query-patterns.md`. `scripts/biomart_query.py` implements the helper.

Done when the table is written, missing mappings and multiplicity are reported, and dataset and release recorded.
