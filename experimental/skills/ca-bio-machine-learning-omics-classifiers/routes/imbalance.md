# Rare positive class

```bash
python scripts/calibration_check.py data.csv --label label --positive X --out calibration_results.tsv
```

Splits train, calibration and test, fits the model unweighted and with `class_weight=balanced`, recalibrates the balanced model on the calibration split, and reports test AUC, Brier, mean predicted over observed prevalence and calibration slope for the three variants. `--model rf` uses a random forest instead of logistic.

- When the probability is the output, change nothing about the training prior: no `class_weight`, `scale_pos_weight` or SMOTE. Choose the operating threshold by cost on a validation fold; the sensitivity resampling buys is recoverable that way.
- Reweighting always moves the probability scale (balanced logistic at 8% prevalence predicted 4.2x the observed rate, Brier 0.177 against 0.070; balanced forest 2.6x). For logistic that is all it changes; for a forest it can change the ranking too, so test it.
- Hard-label problem only: use class weights or SMOTE inside the training fold with an `imblearn` Pipeline (not sklearn's), and recalibrate on held-out data at natural prevalence before reading any risk. `imblearn` was not available when this Skill was recut, so that path is untested here.
- Report AUPRC and MCC under imbalance; accuracy and ROC-AUC mislead when positives are rare. A 95%-accurate model can be useless.
- A calibration split with few events gives a poor calibrator; the script refuses fewer than 3 in the rarer class.

Done when the answer states the prior-shift factor (mean predicted over observed) and the threshold or recalibration used.
