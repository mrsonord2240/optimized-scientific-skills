# Interpret Bayes factors

```bash
python scripts/interpret_bagel.py bayes_factor.txt --screen-type dropout -o calls.tsv
```

Default essential threshold is BF > 6. Only label negative-BF candidates for a screen designed for enrichment; configure assay control genes and calibrate thresholds on this screen's precision-recall. See `references/interpret-calls.md` and `references/method-comparison.md`.

Done when calls, controls and screen-specific threshold are reported.
