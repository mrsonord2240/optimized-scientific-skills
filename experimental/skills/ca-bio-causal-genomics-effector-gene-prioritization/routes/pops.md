# PoPS ranking

```bash
python scripts/pops_run.py --pops_repo /path/to/pops --magma_prefix genomewide_gene --gene_annot_path gene_annot.txt --feature_mat_prefix features --control_features_path control.features --out_prefix pops_out
```

Use genome-wide multi-chromosome MAGMA output, not a single locus; the wrapper defaults to one feature chunk. Set --num_feature_chunks for the actual matrix. Interpret as a locus-relative similarity ranking and inspect feature bias. See `references/pops.md`.

Done when per-gene scores and feature and input scope are reported.
