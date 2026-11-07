# Replicate-aware pseudobulk

```bash
leafcutter_ds.R -i 5 -g 3 -c 20 -o pb_ds pb_counts.txt.gz pb_groups.txt
```

First aggregate junction counts by sample and cell type using `references/pseudobulk.md`; input must be Leafcutter clustered intron counts, not arbitrary junction rows. Require at least three independent samples per group and adequate junction support. rMATS instead needs per-sample-by-cell-type BAMs. Do not count cells as biological replicates.

Done when replicated group results and support thresholds are reported.
