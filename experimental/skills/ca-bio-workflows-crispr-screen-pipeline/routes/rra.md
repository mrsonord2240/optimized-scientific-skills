# Two-condition screen with MAGeCK RRA

```bash
python scripts/rra.py experiment.count.txt essentiality_rra treatment=Day14_r1,Day14_r2,Day14_r3 control=Day0
```

`treatment=` and `control=` are required and must be sample labels from the count table: Day0 or plasmid for dropout, the vehicle arm for a drug screen. Writes `essentiality_rra.gene_summary.txt`, `_negative_hits.csv` (depleted), `_positive_hits.csv` (enriched) and `_volcano.png`.

- A hit needs `fdr=0.05` and `|lfc| >= 0.5`; pass the request's own cutoffs, `lfc=0` to filter on FDR alone.
- Median normalization breaks when more than 40% of guides change (everything significant at FDR 0.01): rerun with `norm=control ctrl_genes=<file>` of non-targeting controls, or use `routes/bagel2.md`.
- RRA is wrong for time courses: `routes/mle.md`.
- Cancer line: the table is the CRISPRcleanR-corrected one.
- Calling a second method for consensus: `routes/consensus.md`.

Done when `gene_summary.txt` exists and the hit tables are read from it.
