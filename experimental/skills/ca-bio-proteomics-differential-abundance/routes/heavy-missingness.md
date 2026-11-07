# Many missing values, proDA

```bash
Rscript scripts/proda_de.R proteinGroups.txt samples.csv Treatment-Control results.csv normalize=median
```

- Missing values enter as left-censored observations. Nothing is imputed.
- proDA needs intensity-dependent missingness (low proteins drop out). If values are missing at random, such as a lost TMT channel, use `routes/protein-table.md` or `routes/psm-counts.md`.
- Proteins with no observed value in a group have blank `diff`, `t_statistic`, `se` and a name in `undetected_in`. List them as "undetected in <group>". Their non-significance is not evidence of no change.
- `batch` in `samples.csv` goes into the model.
- Other options and the input layout are those of `routes/protein-table.md`.

Done when `results.csv` exists and the undetected proteins are listed apart from the significant ones.
