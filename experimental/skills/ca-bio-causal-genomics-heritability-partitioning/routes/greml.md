# Individual-level GREML

```bash
bash scripts/gcta_greml.sh genotype_prefix phenotype.txt result
```

Input phenotype has FID, IID and phenotype. Inspect relatedness, ancestry covariates and sample size; see `references/gcta-greml.md`. The source suggests a PCGC option for rare case-control studies; verify the installed GCTA interface before using it.

Done when h2, SE and model assumptions are reported.
