# Alpha and beta metrics

```bash
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza --p-sampling-depth DEPTH --m-metadata-file metadata.tsv --output-dir core-metrics-results
```

Replace depth with a justified integer and verify the tree covers the retained features. Samples below depth are dropped; report them. The source also provides `examples/core_metrics_qiime2.sh` as a fixed-path demonstration with SEPP and both weighted and unweighted tests; its depth is a placeholder and must be chosen after inspection. `examples/diversity_analysis.R` uses bundled data and is not a user-data CLI. Adapt its alpha, UniFrac, PERMANOVA and dispersion code for actual data and design. PERMANOVA and dispersion tests answer different questions; significant dispersion makes location interpretation ambiguous. Restricted permutations must match the exchangeability structure and the effect being tested; within-subject restrictions cannot test a factor that never varies within subject. See [interpretation](../references/interpretation.md).

Done when metric, tree, depth, exclusions, permutation design and dispersion results are reported.
