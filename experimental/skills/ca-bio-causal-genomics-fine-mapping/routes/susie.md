# Summary-statistic SuSiE

```bash
Rscript scripts/susie_rss.R locus_sumstats.tsv locus_ld.tsv 50000 finemap 10
```

Input columns are SNP, BETA, SE; LD is a signed named square matrix in identical SNP order. Supply actual N, using effective N = 4/(1/Ncase+1/Ncontrol) for case-control statistics. The wrapper computes lambda and kriging diagnostics and stops above 0.10; lambda 0.05-0.10 is marginal, not high-confidence. It writes PIPs and the fit object. Review convergence, credible-set coverage and purity; typical common-variant purity threshold is 0.5. Increase L if all components are used. See `references/diagnostics.md` and `references/reporting.md`. HLA and inversion regions need specialized handling.

Done when diagnostics and credible sets are reviewed and reported with LD source and assumptions.
