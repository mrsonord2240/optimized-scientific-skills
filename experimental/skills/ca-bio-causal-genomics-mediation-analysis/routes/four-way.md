# Four-way decomposition

```bash
Rscript scripts/cmaverse_4way.R data.csv exposure mediator outcome age,sex result.rds
```

This is the source logistic outcome and linear mediator example with interaction, exposure contrast 0 to 1, mediator fixed at 0 and 1000 bootstrap iterations. Confirm these choices and rare-outcome assumptions fit the data; adapt the source method for other scales or exposure levels. See `references/cmaverse-4way.md`. Assess sensitivity separately.

Done when all decomposition terms, assumptions and uncertainty are reported.
