# Multi-signal coloc

```bash
Rscript scripts/coloc_susie.R gwas.tsv eqtl.tsv ld.tsv --gwas-type cc --gwas-s 0.3 --gwas-n 50000 --eqtl-type quant --eqtl-sdy 1 --eqtl-n 500 --L 10 --out coloc_susie_out
```

Tables must be harmonized and share SNP order with the named signed Pearson r matrix. Prefer in-sample or ancestry-matched LD. The wrapper checks each trait's z/LD lambda against 0.05 and gates hg38 MHC and inversion. See `references/coloc-susie.md`. A missing overlapping credible-set pair is not proof of absence.

Done when per-pair posteriors and LD checks or a no-overlap result are reported.
