# Alignment confidence

```bash
muscle -align sequences.fasta -stratified -output ensemble.efa
```

For column-specific conclusions compare independent aligners or score an ensemble and flag unstable columns. `-stratified` and `-diversified` work with `-align`, not `-super5`; see `references/confidence.md` for large-dataset ensembles, TCS, thresholds and validation. Apply tool-specific confidence thresholds rather than transferring one uncalibrated.

Done when requested per-column uncertainty is quantified and unreliable positions are flagged for downstream work.
