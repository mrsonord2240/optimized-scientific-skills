# Several screens in one library (JACKS)

```bash
python <JACKS clone>/jacks/run_JACKS.py experiment.count.txt replicatemap.txt guidemap.txt \
    --rep_hdr Replicate --sample_hdr Sample --ctrl_sample_hdr Control \
    --sgrna_hdr sgRNA --gene_hdr Gene --outprefix jacks_out --apply_w_hp
```

`replicatemap.txt` is tab-separated with columns `Replicate` (a count-table column), `Sample` (the cell line) and `Control` (the `Sample` name of the control pool; the plasmid replicate's row has that name as both `Sample` and `Control`). `guidemap.txt` has columns `sgRNA` and `Gene`, the same names the count table uses for `--sgrna_hdr` and `--gene_hdr`.

- JACKS is a GitHub clone and the script sits in its `jacks/` folder, not the clone root; the PyPI package named `jacks` is unrelated.
- Output: `jacks_out_gene_JACKS_results.txt` (gene effect) and the guide efficacy table.
- Time course: `routes/mle.md`.

Done when `jacks_out_gene_JACKS_results.txt` exists.
