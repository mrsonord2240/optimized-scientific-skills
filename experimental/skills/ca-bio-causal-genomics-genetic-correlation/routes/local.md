# LAVA local correlation

```bash
Rscript scripts/lava_local_rg.R input.info.txt sample.overlap.txt 1kg_EUR_chr blocks.locfile trait1,trait2 lava_out
```

Match build and ancestry across summary statistics, reference and loci. The wrapper requires univariate local h2 p < 0.05/N_loci in both traits before the bivariate test and applies Bonferroni correction. For rho-HESS use `references/hess.md`; LAVA partial correlations require an adapted call.

Done when tested loci, filtered loci, corrected local results and overlap assumptions are reported.
