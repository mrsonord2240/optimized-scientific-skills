# Weighting and coevolution

```bash
python scripts/weights.py alignment.fasta weights.json
```

The source position-based weights skip every column containing a gap; they sum to one and are not an effective sequence count. Neff uses the source identity threshold of 0.62. For gappy alignments or alternative estimators read `references/weighting-neff.md`. For MI-APC, depth guards and shuffled-null interpretation read `references/coevolution-mi-apc.md`; its implementation is `examples/mi_apc.py`.

To calculate Neff independently when no gap-free columns exist, run `python examples/neff.py alignment.fasta` and record its estimator and threshold.

Done when the selected estimator, its assumptions, and resulting weights or scores are reported.
