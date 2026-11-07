# PSM or peptide counts available, DEqMS

```bash
Rscript scripts/deqms_de.R proteinGroups.txt samples.csv Treatment-Control results.csv count_col=Peptides normalize=median
```

- Counts from a separate file instead: `count=counts.csv` with columns `protein`, `count`.
- Count: PSMs for TMT, peptides for label-free. Multi-batch TMT: the minimum across batches, one value per protein.
- Every tested protein needs a count of at least 1; the script stops otherwise.
- Read `sca.adj.pval` and `sca.P.Value`, not the limma columns.
- Other options and the input layout are those of `routes/protein-table.md`.

Done when `results.csv` exists and you have answered from `sca.adj.pval`.
