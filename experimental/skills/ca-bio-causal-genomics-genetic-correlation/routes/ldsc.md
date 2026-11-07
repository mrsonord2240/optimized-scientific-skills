# Cross-trait LDSC

```bash
bash scripts/ldsc_crosstrait_rg.sh trait1.tsv trait2.tsv 250000 180000 eur_w_ld_chr/ w_hm3.snplist rg_out
```

Supply actual effective sample sizes and ancestry-matched LD scores. The wrapper munges both inputs and warns when mean chi-square is <= 1.02. Report rg, SE, p and cross-trait intercept; the intercept can absorb sample overlap under LDSC assumptions. The wrapper's parsed trigger is heuristic; inspect the log rather than treating it as a causal conclusion.

Done when the log and uncertainty, ancestry and overlap context are reported.
