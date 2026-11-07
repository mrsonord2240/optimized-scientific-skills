# regenie set tests

```bash
bash examples/rare_variant_test.sh array_prefix exome_prefix pheno.txt covar.txt annot.txt sets.txt masks.txt output_dir
```

The source script checks mask concordance, fits a binary-trait step-1 LOCO predictor using QC'd common variants, and runs step 2 on the rare-variant set with AAF bins 0.001 and 0.01, `skato,acato`, max mask construction and approximate Firth at pThresh 0.05. Supply PLINK prefixes, phenotype and covariate files, annotation, set-list and mask-definition files. It creates the specified directory if absent and uses fixed output names, so use a fresh directory. Annotation rows are variant, gene, annotation; set rows are gene, chromosome, position, comma-separated variant IDs; mask rows are mask name and comma-separated annotations. The source notes that `.all` is implicit and singleton masks require a separate option. Check installed help: the source example targets regenie 4.1.3 and interfaces were not verified here.

**Done when:** the mask concordance report, `gene_tests_*.regenie` outputs, masks, null provenance and per-gene/per-mask testing family are retained and reviewed.
