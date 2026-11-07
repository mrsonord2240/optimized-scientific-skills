# Interpret corrected labeling

```bash
python scripts/correct_csv.py intensities.csv FORMULA 0.99 corrected.csv
```

This obtains a corrected MID for one input vector; see [correction](correction.md) for required metadata. Mean enrichment is the atom-weighted fraction, but neither it nor a single MID uniquely determines flux. The source suggests a descriptive plateau screen requiring at least four timepoints and the last three consecutive enrichment changes below 0.02; this is a heuristic, not proof of isotopic or metabolic steady state. For transients, uncertainty and flux inference, use an appropriate model and independent validation; no time-course or flux-fitting CLI is supplied here. Pool concentration and labeling answer different questions.

Done when timepoints, correction, uncertainty, steady-state assumptions and limits of inference are stated.
