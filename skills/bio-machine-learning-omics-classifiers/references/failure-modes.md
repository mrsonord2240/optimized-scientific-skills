# Omics Classifier Failure Modes

Classifier failure modes: batch-artifact AUC, uncalibrated RF/boosting probabilities, class-weight/SMOTE misuse, tuning data reported as test, OOB-as-test.

## Per-Method Failure Modes

### Beautiful AUC that is a batch artifact
- **Trigger:** Cases and controls processed in different batches/sites/times.
- **Mechanism:** The classifier exploits the cleaner technical signal; even nested CV is optimistic, and ComBat cannot rescue a confounded design (Soneson 2014).
- **Symptom:** Near-perfect CV AUC; collapse on an independent cohort; the model predicts batch easily.
- **Fix:** Detect with the batch-prediction check; use leave-one-batch-out (`scripts/batch_checks.py`, per-batch AUC); fix at design (balance batches across outcome).

### RF/boosting probabilities trusted as risks
- **Trigger:** Reading `predict_proba` from RF or XGBoost as a calibrated risk.
- **Mechanism:** RF bounds votes away from 0 and 1; classic boosting is sigmoid-distorted (Niculescu-Mizil 2005) and log-loss GBDT can overfit to overconfident extremes.
- **Symptom:** Good AUC, reliability curve far from diagonal.
- **Fix:** Recalibrate on a disjoint fold; or prefer logistic when the probability matters.

### Class weights / SMOTE for a risk model, or SMOTE before the split
- **Trigger:** `class_weight='balanced'`, `scale_pos_weight` or resampling to "fix" imbalance for a probability model, or resampling before the CV split.
- **Mechanism:** Reweighting and resampling change the training prior, so predicted minority risk inflates (4.2x the prevalence for balanced logistic at 8%, see SKILL.md); synthetic points derived from test samples leak when resampled before the split (van den Goorbergh 2022; Carriero 2025).
- **Symptom:** Risks systematically too high (AUC unchanged for logistic, can shift for RF: 0.776 to 0.812 in SKILL.md); inflated CV performance if resampled before the split.
- **Fix:** Do not reweight or resample risk models; tune the threshold. For hard-label problems use an `imblearn` Pipeline (train-fold only) and recalibrate on held-out data.

### Early-stopping or calibration data reported as the result
- **Trigger:** Reporting the XGBoost early-stopping validation metric, or a Brier computed on the calibration set.
- **Mechanism:** Both sets were used to choose the model (rounds, calibrator); a small set saturates, and isotonic on few samples returns only 0 and 1.
- **Symptom:** Validation AUCPR near 1 with a much lower test AUC; Brier 0.000; three or fewer distinct probabilities.
- **Fix:** Report on a third, untouched split; thresholds and sample sizes are in SKILL.md.

### OOB error read as an unbiased test estimate
- **Trigger:** Reporting RF out-of-bag error after selecting features on the full data.
- **Mechanism:** Selection leaked; OOB then reflects the contaminated feature set.
- **Symptom:** Optimistic OOB; external collapse.
- **Fix:** Selection inside CV; estimate performance by nested CV.
