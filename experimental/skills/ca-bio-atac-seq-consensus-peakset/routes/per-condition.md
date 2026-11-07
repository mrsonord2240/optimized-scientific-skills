# Consensus per condition, then union

```bash
bash scripts/per_condition_consensus.sh consensus_out 2 \
  cond1=cond1_rep1.narrowPeak,cond1_rep2.narrowPeak,cond1_rep3.narrowPeak \
  cond2=cond2_rep1.narrowPeak,cond2_rep2.narrowPeak,cond2_rep3.narrowPeak
```

Arguments: output directory, the minimum replicates per condition (peak present in at least that many), then one `name=file1,file2,...` per condition. It writes `<cond>_consensus.bed` per condition and `all_conditions_consensus.bed`.

- The union is variable width. For differential counting, re-center the inputs to fixed width first, or run `routes/iterative-overlap.md` on the pooled replicates instead.
- A union much larger than any single condition's consensus is the expected result of real condition-specific peaks.

Done when `all_conditions_consensus.bed` exists; next `routes/count.md`.
