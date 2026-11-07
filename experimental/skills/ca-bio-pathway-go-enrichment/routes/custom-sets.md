# Your own gene sets

```bash
Rscript scripts/ora_go.R genes.txt universe.txt results.csv term2gene=sets.tsv keytype=SYMBOL
```

`sets.tsv`: two columns, term then gene, no header, in the same ID type as `genes.txt` (here SYMBOL). The script runs `enricher`, the same hypergeometric test as `enrichGO`, with the same explicit universe. `ont` and `simplify` do not apply.

- With a DE table use `de=de_results.csv` in place of the two gene files (see `routes/de-table.md`).
- Sets with fewer than 10 or more than 500 genes inside the universe are dropped (`minGSSize`, `maxGSSize`).

Done when `results.csv` exists and the answer reports `p.adjust`, `Count` and the universe size.
