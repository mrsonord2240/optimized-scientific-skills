# scQuint plate junctions

```bash
python scripts/scquint_groups.py Solo.out/SJ/raw annotation.gtf.gz cell_labels.tsv neuron glia scquint
```

The label TSV uses cell ID as its first column and cell_type as a column; select the user's actual groups. The wrapper follows `references/scquint.md` with 10-cell thresholds and writes group and intron tables. The source version expects chr-prefixed junctions and an unprefixed GTF; inspect annotation and report the threshold choice. Do not use 10X 3-prime or 5-prime for this route. Group comparisons need donor-aware interpretation.

Done when supported intron-group results are produced and replication limitations are reported.
