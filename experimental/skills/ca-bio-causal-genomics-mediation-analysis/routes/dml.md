# Double-ML mediation

```bash
Rscript scripts/meddml.R data.csv outcome treatment mediator age,sex result.rds
```

Use numeric named covariates and verify overlap, treatment and mediator definitions. The source example uses trim=0.05 and order=1; inspect trimmed counts and assumptions in `references/time-varying-and-dml.md`.

Done when direct and indirect effects, uncertainty and trimming are reported.
