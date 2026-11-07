# Combine predictions

```bash
python scripts/combine_predictions.py --abc abc_threshold.tsv --re2g re2g_threshold.tsv.gz --loops loops.bedpe --out combined.tsv
```

Both tables must use the same candidate regions and each method's calibrated threshold; omit --loops if unavailable. Intersections are higher-priority hypotheses and single-method calls exploratory. Compare to independent validation data where available; loop convention is FDR < 0.05 and at least five contacts.

Done when overlap and any loop support are reported with thresholds and limitations.
