# Gene and ortholog metadata

```bash
datasets summary gene symbol BRCA1 --ortholog Mammalia --as-json-lines | dataformat tsv gene --fields gene-id,symbol,tax-name,description,chromosomes > brca1_mammals.tsv
```

Replace the symbol and clade. `--ortholog` requires an explicit value such as `all`; `--taxon` resolves a gene at species level and does not accept a clade for this purpose. The source describes NCBI's representative-per-species set as limited to vertebrates and insects; it is not a full co-ortholog analysis. For single species use `datasets summary gene symbol BRCA1 --taxon human --as-json-lines`. Check field names against installed help. See [gene details](../references/gene-orthologs.md).

Done when the TSV contains the intended genes and species; inspect count and identity rather than exit status alone.
