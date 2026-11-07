# Ortholog table

```bash
python scripts/query.py --dataset hsapiens_gene_ensembl --attributes ensembl_gene_id,external_gene_name,mmusculus_homolog_ensembl_gene,mmusculus_homolog_orthology_type,drerio_homolog_ensembl_gene,drerio_homolog_orthology_type --filter chromosome_name --value 17 --output orthologs.tsv
```

Replace chromosome and species attributes as requested. Column display names and orthology classifications can differ by release; inspect before filtering one-to-one relationships. For version-pinned R and other examples read `references/query-patterns.md`.

Done when the requested ortholog rows and their type, dataset and release are reported.
