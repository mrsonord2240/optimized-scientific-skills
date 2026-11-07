# Calibration

```bash
python scripts/calibration.py predictions.csv --label label --prob prob
```

`predictions.csv`: one row per held-out sample with the true label and the predicted probability, from samples the model was not trained on. The script prints the reliability curve (equal-mass bins, `--bins 10`), Brier score against the prevalence-only Brier, and the calibration intercept and slope from a logistic fit of the label on the logit of the probability.

- Slope below 1 means overconfident predictions (overfitting); intercept away from 0 means predicted risk is systematically high or low.
- Do not report one Expected Calibration Error as the verdict. Equal-width ECE is biased.
- AUC cannot show miscalibration; report the curve, slope and Brier alongside it.

Recalibrate a fitted model only on a fold disjoint from both training and evaluation data. Platt (`sigmoid`) for small calibration sets, isotonic for hundreds of points or more:

```python
from sklearn.calibration import CalibratedClassifierCV
from sklearn.frozen import FrozenEstimator   # scikit-learn >= 1.6; cv='prefit' no longer exists

calibrated = CalibratedClassifierCV(FrozenEstimator(fitted_model), method='isotonic')
calibrated.fit(X_cal, y_cal)                  # X_cal disjoint from train and test
p_test = calibrated.predict_proba(X_test)[:, 1]
```

Recalibrating on the evaluation data makes the curve look perfect and is leakage.

Done when you have reported slope, intercept, Brier and the reliability curve from data the model and the recalibrator never saw.
