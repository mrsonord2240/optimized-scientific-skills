# Time course, several conditions, or batches (MAGeCK MLE)

```bash
mageck mle --count-table experiment.count.txt --design-matrix design.txt \
    --output-prefix timecourse_mle --norm-method median --permutation-round 10
```

`design.txt` is tab-separated: first column `Samples`, then a `baseline` column of 1s, then one 0/1 column per condition or covariate. Samples are the count table's labels.

- Batch: add one 0/1 column per extra batch to the design matrix. Never correct batch with ComBat on counts first.
- Permutation rounds default to 2: use at least 10, then recheck borderline `fdr` calls near 0.05 against `wald-fdr`. Do not trust a default-round boundary call.
- Genome-scale runs take hours.
- Results are per-condition beta scores in `timecourse_mle.gene_summary.txt`.

Done when `gene_summary.txt` exists and borderline genes were rechecked.
