# High-dimensional mediators

```bash
Rscript scripts/hima_ewas.R pheno.csv mediators.csv "outcome ~ exposure + age + sex" out.csv gaussian DBlasso 0.05 8
```

Rows must align between phenotype and mediator files. Set mediator type for the actual distribution and use HIMA >=2.3.0. The script handles phenotype missingness and categorical covariates; inspect mediator completeness and sample alignment. See `references/hima-ewas.md`.

Done when screened mediators, FDR, missingness and assumptions are reported.
