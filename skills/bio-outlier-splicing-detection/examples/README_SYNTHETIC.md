# Synthetic audit cohort

`01_make_synth_cohort.py` makes a deterministic, entirely synthetic 30-sample RNA-seq-like cohort for exercising the FRASER, OUTRIDER, and LeafcutterMD preparation paths. It requires Python with `pysam` and writes BAMs, BAM indexes, `gene_counts.tsv`, `planted_truth.tsv`, and `genes.tsv`.

```bash
python examples/01_make_synth_cohort.py /tmp/outlier-splicing-synthetic
```

Expected truth: four aberrant-splicing events (S05, S12, S20, S25), two expression events (S15, S28), and one deliberately global tissue-mismatch sample (S29). It is test data only; do not use it for biological or clinical inference.
