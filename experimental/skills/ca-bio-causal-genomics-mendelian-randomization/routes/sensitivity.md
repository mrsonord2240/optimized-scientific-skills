# MR-PRESSO and SIMEX

```bash
Rscript scripts/mr_presso_outliers.R --dat mr_out/harmonised.tsv --nb 10000 --seed 42 --signif 0.05 --out mr_out/presso_outliers.txt
```

Inspect global, outlier and distortion results; outlier p-values are already adjusted. Choose enough simulations for the number of instruments. For Egger NOME violation use `Rscript scripts/simex_egger.R --dat mr_out/harmonised.tsv --B 1000 --seed 1`; it prints slopes. See `references/mr-presso.md` and `references/simex-egger-nome.md`. Outlier tests do not rule out correlated pleiotropy.

Done when requested diagnostics, Monte Carlo precision and assumptions are reported.
