# Annotated two-group BAMs

```bash
bash scripts/diff_splicing_rmats.sh
```

Edit the configuration at the top first: GTF, one-line comma-separated BAM lists, strandedness, output and read settings. The script detects read type and length, uses cstat 0.01, and checks for nonempty SE output with FDR. For paired samples, PAIRADISE and corresponding ordered pairs are required; see `references/rmats.md`. Default reporting FDR < 0.05, |delta PSI| > 0.10 with >=10 junction reads in at least half the replicates of each group. rMATS sign is group 1 minus group 2.

Done when the event tables and requested filtered results are interpreted.
