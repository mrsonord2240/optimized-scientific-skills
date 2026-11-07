# Observational mediation

```bash
Rscript scripts/observational_mediation.R data.csv exposure mediator outcome age,sex result.rds probit
```

Supply actual column names and confounders, or - for none. The wrapper fits a linear mediator and linear or probit outcome, bootstrap mediate with 1000 simulations and rho sensitivity. For final precision the source recommends 5000 simulations; modify deliberately and document. Probit is used for binary outcomes because medsens does not support logit. Check causal assumptions and missing-data handling; see `references/sensitivity-and-reporting.md`.

Done when direct, indirect and total effects, intervals and sensitivity are reported.
