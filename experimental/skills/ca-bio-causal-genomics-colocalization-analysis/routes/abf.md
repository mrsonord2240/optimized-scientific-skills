# Single-causal coloc

```bash
Rscript scripts/coloc_abf.R gwas.tsv eqtl.tsv --gwas-type cc --gwas-s 0.30 --gwas-n 50000 --eqtl-type quant --eqtl-sdy 1 --eqtl-n 500 --p12 5e-6 --out coloc_out
```

First harmonize with `Rscript scripts/harmonise.R gwas_raw.tsv eqtl_raw.tsv harmonised.tsv` and prepare separate aligned tables in identical SNP order, as in `references/allele-harmonisation.md`. Supply real sample sizes, case fraction and trait type; sdY=1 assumes standardized quantitative trait. The wrapper defaults p1=p2=1e-4, p12=5e-6, and gates hg38 MHC and inversion. It writes sensitivity and posterior tables. For p-value/MAF-only input use coloc directly after checking its installed interface. Use lower p12 for trans and report its sensitivity; see `references/priors-and-thresholds.md`. Single-causal assumptions can mislead when multiple signals exist.

Done when PP.H0-H4, p12 sensitivity and limitations are reported.
