# TOBIAS differential footprints

```bash
bash scripts/run_tobias.sh cond1.bam cond2.bam consensus.bed genome.fa blacklist.bed motifs.pfm tobias_out 16
```

Use deduplicated, filtered, mitochondrial-read-removed BAMs and identical peaks and blacklist for both conditions. Source guidance treats fewer than roughly 50M nuclear reads as low power. The script runs bias correction, ScoreBigwig and BINDetect; CTCF controls compare all, bound and unbound profiles. It exits 3 without CTCF unless `ALLOW_NO_CTCF=1`, and exits 4 if the source bias-correlation check exceeds `BIAS_R_MAX` (default 0.2). Use another positive control if CTCF is absent. `PMAX` defaults to 0.05. See `references/interpretation.md` for direction, per-TF limitations and motif choices; see `references/inputs.md` for NFR filtering and single-condition commands.

Done when the calls, differential table and controls are inspected and tool, genome, motif, threshold and bias model are reported.
