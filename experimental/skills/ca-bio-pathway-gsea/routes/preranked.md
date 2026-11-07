# Preranked analysis

```bash
Rscript scripts/gsea_table.R ranking.csv term2gene.csv result.csv result.rds
```

Ranking CSV requires unique nonmissing `gene,stat` with a signed finite statistic for the full tested universe. TERM2GENE CSV requires `term,gene` in the same identifier namespace. The wrapper sorts descending and uses source settings exponent 1, set sizes 10–500, `eps=0`, BH, p-value cutoff 0.05 and seed 123; it writes table and RDS and may overwrite outputs. A raw p-value is unsigned; identify the contrast direction and document any statistic transformation, ties and ID exclusions. Preranked gene-level nulls may not preserve sample-level gene correlations; interpret as screening and compare an appropriate matrix-based method where possible.

Done when mapping coverage, ranking, contrast, collection version, settings, adjusted results and leading edges are reviewed.
